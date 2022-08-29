from importlib import resources
from pathlib import Path
import pandas as pd

from pyomo.environ import (
    Var,
    NonNegativeReals,
    Binary,
    Param,
    RangeSet,
    Constraint,
    Expression,
)

from idaes.core.base.process_base import declare_process_block_class
from idaes.models.unit_models import SkeletonUnitModelData
from pyomo.common.config import ConfigValue

# TODO: Contact Sandeep for reference. 
# HHV of NG = 22499.17034 btu/lb = 0.0496 MMBtu/kg
NG_HHV = 0.0496
HR_TO_SEC = 3600


def get_lmp_data(m, 
                 dataset="NREL",
                 location="CAISO",
                 carbon_tax=100,
                 princeton_case="BaseCaseTax"):
    """
    This function reads and appends LMP data to the model.
    """
    with resources.path(
        "idaes.power_generation.flowsheets.direct_fired_cycle", 
        "FLECCS_Price_Series_Data_01_20_2021.xlsx"
    ) as p:
        path_to_file = Path(p).resolve()

    if dataset == "NREL":
        raw_data = pd.read_excel(
            path_to_file,
            sheet_name="2035 - NREL",
        )
        column_name = 'MiNg_$' + str(carbon_tax) + '_' + location
        # We will take the price signal for day 365 same as that for day 364
        num_days = 365

    elif dataset == "Princeton":
        raw_data = pd.read_excel(
            path_to_file,
            sheet_name="2030 - Princeton")
        column_name = princeton_case
        num_days = 365

    m.set_time = RangeSet(num_days * 24)
    price_all = raw_data[column_name].tolist()

    # Set prices lower than $0.001/MWh to zero to avoid numerical issues
    price_all = [i if i > 1e-3 else 0 for i in price_all]

    # LMP is divided by 1000 to formulate the objective in thousands of dollars.
    # This helps in keeping the scale of the CAPEX in O(1e6), otherwise it would be O(1e9).
    m.LMP = Param(
        m.set_time,
        initialize={t + 1: lmp / 1000 for (t, lmp) in enumerate(price_all)},
        doc="Locational Marginal Prices [in 1000$/MWh]"
    )


@declare_process_block_class("DFCDesign")
class DFCDesignData(SkeletonUnitModelData):
    """
    This class contains variables/parameters related to the design of the direct fired cycle
    """
    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare("capacity_range", ConfigValue(
        default=(500, 850),
        doc="Range for the capacity of the DFC [in MW]",
    ))

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        self.capacity = Var(
            within=NonNegativeReals,
            initialize=838.11322145,
            doc="Capacity of the power plant [in MW]",
        )

        self.ng_flow = Var(
            within=NonNegativeReals,
            initialize=28.87,
            doc="Natural gas flowrate at full load [in kg/s]",
        )

        # Define a variable that informs whether the plant needs to built or not. 
        self.build_dfc = Var(
            within=Binary,
            doc="1: Plant is built, 0: Plant is not built",
        )

        # Bound the capcity of the plant in the specified range
        capacity_range = self.config.capacity_range
        self.capacity_lb_con = Constraint(expr=self.capacity >= self.build_dfc * capacity_range[0])
        self.capacity_ub_con = Constraint(expr=self.capacity <= self.build_dfc * capacity_range[1])

        # Compute the natural gas flowrate required at maximum capacity. 
        # FIXME: Assuming a linear relation for now. Update the equation when we have more data
        self.ng_flow_requirement = Constraint(
            expr=self.ng_flow == self.capacity * (28.87 / 838.11322145),
            doc="Computes the natural flowrate required [in kg/s] at full load",
        )

        self.o2_flow = Expression(
            expr=3.785632948 * self.ng_flow,
            doc="Flowrate of oxygen [in kg/s] required at full load",
        )

        # Assumuing that the capex varies linearly with capacity
        self.capex = Expression(
            expr=1128855 * (self.capacity / 838.11322145),
            doc="CAPEX of the power cycle [in 1000$]",
        )

        # Calculating the FOM as 3.157% of the CAPEX. The percentage value is obtained from
        # 35,641.27 (FOM) / 1,128,855 (CAPEX). The FOM is also assumed to vary linearly with the capacity.
        self.fom = Expression(
            expr=0.03157 * self.capex,
            doc="Fixed O&M cost [in 1000$/year]",
        )


@declare_process_block_class("DFCOperation")
class DFCOperationData(SkeletonUnitModelData):
    """
    This class contains the model for the operation of the power cycle
    """

    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare("design_blk", ConfigValue(
        doc="Pointer to the object containing the DFC design information",
    ))
    CONFIG.declare("operating_range", ConfigValue(
        default=(0.2, 1),
        doc="Off-design operating range. Default: (0.2 * full load, full load)",
    ))

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        design_blk = self.config.design_blk

        # Declare variables
        self.power = Var(
            within=NonNegativeReals,
            doc="Total power generated [in MW]",
        )
        self.ng_flow = Var(
            within=NonNegativeReals,
            doc="Natural gas flowrate [in kg/s]",
        )
        self.norm_power = Var(
            within=NonNegativeReals,
            bounds=(0, 1),
            doc="Power normalized with the design value [-]",
        )
        self.norm_ng_flow = Var(
            within=NonNegativeReals,
            bounds=(0, 1),
            doc="Natural gas flowrate normalized with the design value [-]",
        )
        self.op_mode = Var(
            within=Binary,
            doc="1: In Operation, 0: Shutdown",
        )
        self.startup = Var(
            within=Binary,
            doc="1: Plant is startup at this hour, 0: Otherwise",
        )
        self.shutdown = Var(
            within=Binary,
            doc="1: Plant is shutdown at this hour, 0: Otherwise",
        )

        # Definition of normalized power and natural gas flow variables
        self.norm_power_definition = Constraint(
            expr=self.norm_power * design_blk.capacity == self.power,
            doc="Definition of the noramlized power variable",
        )
        self.norm_ng_flow_definition = Constraint(
            expr=self.norm_ng_flow * design_blk.ng_flow == self.ng_flow,
            doc="Definition of the normalized natural gas flowrate variable",
        )

        # Power production as a function of natural gas flowrate
        # We construct a surrogate model of the form (1 - norm_power) = m * (1 - norm_ng_flow).
        # This way, the natural gas requirement is exact at full load.
        # Rearranging the equation yields norm_power = m * norm_ng_flow + 1-m 
        self.power_production = Constraint(
            expr=self.norm_power == 1.2845 * self.norm_ng_flow - 0.2845 * self.op_mode,
        )

        # Ensure that power production is within P_min and P_max
        operating_range = self.config.operating_range
        self.power_production_lb = Constraint(
            expr=operating_range[0] * self.op_mode <= self.norm_power,
        )

        self.power_production_ub = Constraint(
            expr=self.norm_power <= operating_range[1] * self.op_mode,
        )

        # Declare a variable to track NG requirement for startup and shutdown
        self.su_sd_ng_flow = Var(
            within=NonNegativeReals,
            doc="Fuel requirement for startup/shutdown [in kg/s]",
        )

        self.total_ng_flow = Expression(
            expr=self.ng_flow + self.su_sd_ng_flow,
            doc="Total NG flow rate requirement [in kg/s]",
        )

        self.o2_flow = Expression(
            expr=3.785632948 * self.total_ng_flow,
            doc="Flowrate of oxygen required [in kg/s]",
        )


@declare_process_block_class("MonoASUDesign")
class MonoASUDesignData(SkeletonUnitModelData):
    """
    This class builds the model for the design of the ASU. Here, we treat the
    system of ASUs as a monolithic ASU.
    """

    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare("o2_flow_range", ConfigValue(
        default=(80, 130),
        doc="Range for the size of the ASU in terms of O2 flow rate [in kg/s]",
    ))

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        self.max_o2_flow = Var(
            within=NonNegativeReals,
            initialize=109.2912232,
            doc="Maximum flowrate of O2 the ASU can produce [in kg/s]",
        )

        self.max_power = Var(
            within=NonNegativeReals,
            initialize=162.1878617,
            doc="Power requirement at maximum capacity [in MW]",
        )

        self.build_asu = Var(
            within=Binary,
            doc="1: ASU is built, 0: ASU is not built",
        )

        # Bound the capcity of the plant in the specified range
        o2_flow_range = self.config.o2_flow_range
        self.o2_flow_lb_con = Constraint(expr=self.max_o2_flow >= self.build_asu * o2_flow_range[0])
        self.o2_flow_ub_con = Constraint(expr=self.max_o2_flow <= self.build_asu * o2_flow_range[1])

        # Relation between the flowrate and the power requirement
        self.power_requirement = Constraint(
            expr=self.max_power == self.max_o2_flow * (162.1878617 / 109.2912232)
        )

        # Assuming that the capex of the ASU varies linearly with size
        self.capex = Expression(
            expr=545522 * (self.max_o2_flow / 109.2912232),
            doc="CAPEX of the ASU unit [in 1000$]",
        )

        # Calculating the FOM as 3.157% of the CAPEX. The percentage value is obtained from
        # 17223.73 (FOM) / 545,522 (CAPEX). The FOM is also assumed to vary linearly with the capacity. 
        self.fom = Expression(
            expr=0.03157 * self.capex,
            doc="Fixed O&M cost [in $1000]",
        )


@declare_process_block_class("MonoASUOperation")
class MonoASUOperationData(SkeletonUnitModelData):
    """
    This class contains the model for the operation of the monolithic ASU
    """

    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare("design_blk", ConfigValue(
        doc="Pointer to the object containing the ASU design information",
    ))
    CONFIG.declare("operating_range", ConfigValue(
        default=(0.3, 1),
        doc="Off-design operating range. Default: (0.3 * full load, full load)",
    ))

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        design_blk = self.config.design_blk

        # Declare variables
        self.power = Var(
            within=NonNegativeReals,
            doc="Power required to produce oxygen [in MW]",
        )
        self.o2_flow = Var(
            within=NonNegativeReals,
            doc="Flowrate of oxygen produced [in kg/s]",
        )
        self.norm_power = Var(
            within=NonNegativeReals,
            bounds=(0, 1),
            doc="Power normalized with the design value [-]",
        )
        self.norm_o2_flow = Var(
            within=NonNegativeReals,
            bounds=(0, 1),
            doc="oxygen flowrate normalized with the design value [-]",
        )
        self.op_mode = Var(
            within=Binary,
            doc="1: In Operation, 0: Shutdown",
        )
        self.startup = Var(
            within=Binary,
            doc="1: ASU is startup at this hour, 0: Otherwise",
        )
        self.shutdown = Var(
            within=Binary,
            doc="1: ASU is shutdown at this hour, 0: Otherwise",
        )

        # Definition of normalized power and natural gas flow variables
        self.norm_power_definition = Constraint(
            expr=self.norm_power * design_blk.max_power == self.power,
            doc="Definition of the noramlized power variable",
        )
        self.norm_o2_flow_definition = Constraint(
            expr=self.norm_o2_flow * design_blk.max_o2_flow == self.o2_flow,
            doc="Definition of the normalized oxygen flowrate variable",
        )

        # Power production as a function of natural gas flowrate
        # We construct a surrogate model of the form (1 - norm_power) = m * (1 - norm_o2_flow).
        # This way, the power requirement is exact at full capacity.
        # Rearranging the equation yields norm_power = m * norm_o2_flow + 1-m 
        self.power_requirement = Constraint(
            expr=self.norm_power == 0.9625 * self.norm_o2_flow + 0.0375 * self.op_mode, 
        )

        # Ensure that the normalized oxygen flowrate is within the admissible operating range
        operating_range = self.config.operating_range
        self.o2_flow_lb = Constraint(
            expr=operating_range[0] * self.op_mode <= self.norm_o2_flow,
        )

        self.o2_flow_ub = Constraint(
            expr=self.norm_o2_flow <= operating_range[1] * self.op_mode,
        )

        self.su_sd_power = Var(
            within=NonNegativeReals,
            doc="Power requirement during startup and shutdown [in MW]",
        )

        self.total_power = Expression(
            expr=self.power + self.su_sd_power,
            doc="Total power required by the ASU [in MW]",
        )


@declare_process_block_class("NLUDesign")
class NLUDesignData(SkeletonUnitModelData):
    """
    This class builds the model for the design of the NLU.
    """

    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare("o2_flow_range", ConfigValue(
        default=(50, 130),
        doc="Range for the size of the NLU in terms of O2 flow rate [in kg/s]",
    ))

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        self.max_o2_flow = Var(
            within=NonNegativeReals,
            initialize=109.291223201103,
            doc="Maximum flowrate of O2 the NLU can liquefy [in kg/s]",
        )

        self.max_power = Var(
            within=NonNegativeReals,
            initialize=143.8024,
            doc="Power requirement at maximum capacity [in MW]",
        )

        self.build_nlu = Var(
            within=Binary,
            doc="1: NLU is built, 0: NLU is not built",
        )

        # Bound the capcity of the plant in the specified range
        o2_flow_range = self.config.o2_flow_range
        self.o2_flow_lb_con = Constraint(expr=self.max_o2_flow >= self.build_nlu * o2_flow_range[0])
        self.o2_flow_ub_con = Constraint(expr=self.max_o2_flow <= self.build_nlu * o2_flow_range[1])

        # Relation between the flowrate and the power requirement
        self.power_requirement = Constraint(
            expr=self.max_power == self.max_o2_flow * (143.8024 / 109.291223201103)
        )

        # Assuming that the capex of the NLU varies linearly with size
        # Obtained from $42,797.4 * 4 = $171,189.6
        self.capex = Expression(
            expr=171189.6 * (self.max_o2_flow / 109.291223201103),
            doc="CAPEX of the ASU unit [in 1000$]",
        )

        # Assuming that the FOM of NLU is 3.157% of CAPEX. This number is calculated as
        # 41,715 (FOM) / 171,189 (CAPEX) 
        self.fom = Expression(
            expr=0.03157 * self.capex,
            doc="Fixed O&M of NLU [in $1000]",
        )


@declare_process_block_class("NLUOperation")
class NLUOperationData(SkeletonUnitModelData):
    """
    This class contains the model for the operation of the NLU
    """
    # NOTE: The NLU model can be simplified by not defining the design variables.
    # But we modeled it in this manner to keep it consistent with the rest of the models.
    # If the solver finds it difficult to solve, then we'll switch to the simpler model.

    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare("design_blk", ConfigValue(
        doc="Pointer to the object containing the NLU design information",
    ))
    CONFIG.declare("operating_range", ConfigValue(
        default=(0.3, 1),
        doc="Off-design operating range. Default: (0.3 * full load, full load)",
    ))

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        design_blk = self.config.design_blk

        # Declare variables
        self.power = Var(
            within=NonNegativeReals,
            doc="Power required to liquify oxygen [in MW]",
        )
        self.o2_flow = Var(
            within=NonNegativeReals,
            doc="Flowrate of liquified oxygen [in kg/s]",
        )
        self.norm_power = Var(
            within=NonNegativeReals,
            bounds=(0, 1),
            doc="Power normalized with the design value [-]",
        )
        self.norm_o2_flow = Var(
            within=NonNegativeReals,
            bounds=(0, 1),
            doc="oxygen flowrate normalized with the design value [-]",
        )

        # Note: Assuming that the NLU does not have any startup/shutdown constraints.
        # If it does, we need to define startup and shutdown binary variables like before.
        self.op_mode = Var(
            within=Binary,
            doc="1: In Operation, 0: Shutdown",
        )

        # Definition of normalized power and natural gas flow variables
        self.norm_power_definition = Constraint(
            expr=self.norm_power * design_blk.max_power == self.power,
            doc="Definition of the noramlized power variable",
        )
        self.norm_ng_flow_definition = Constraint(
            expr=self.norm_o2_flow * design_blk.max_o2_flow == self.o2_flow,
            doc="Definition of the normalized oxygen flowrate variable",
        )

        # Power production as a function of natural gas flowrate
        self.power_requirement = Constraint(
            expr=self.norm_power == self.norm_o2_flow,
        )

        # Ensure that the normalized oxygen flowrate is within the operating_range
        operating_range = self.config.operating_range
        self.o2_flow_lb = Constraint(
            expr=operating_range[0] * self.op_mode <= self.norm_o2_flow,
        )

        self.o2_flow_ub = Constraint(
            expr=self.norm_o2_flow <= operating_range[1] * self.op_mode,
        )


@declare_process_block_class("OxygenTankDesign")
class OxygenTankDesignData(SkeletonUnitModelData):
    """
    This class builds the model for the design of the tank.
    """

    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare("tank_size_range", ConfigValue(
        default=(2000, 400000),
        doc="Range for the tank capacity [in tons]",
    ))

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        self.tank_capacity = Var(
            within=NonNegativeReals,
            initialize=4000,
            doc="Maximum amount of oxygen the tank can hold [in tons]",
        )

        self.build_tank = Var(
            within=Binary,
            doc="1: The storage tank is built, 0: Otherwise",
        )

        # Bound the capcity of the tank in the specified range
        tank_size_range = self.config.tank_size_range
        self.tank_capacity_lb_con = Constraint(expr=self.tank_capacity >= self.build_tank * tank_size_range[0])
        self.tank_capacity_ub_con = Constraint(expr=self.tank_capacity <= self.build_tank * tank_size_range[1])

        self.capex = Expression(
            expr=0.98167214 * self.tank_capacity + 2779.90543786 * self.build_tank,
            doc="CAPEX of the oxygen storage tank [in 1000$]",
        )

        self.fom = Expression(
            expr=0.0317 * self.capex,
            doc="FOM of the tank [in $1000]",
        )


@declare_process_block_class("OxygenTankOperation")
class OxygenTankOperationData(SkeletonUnitModelData):
    """
    This class contains the model for the operation of the oxygen tank
    """

    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare("design_blk", ConfigValue(
        doc="Pointer to the object containing the tank design information",
    ))
    CONFIG.declare("o2_flow_range", ConfigValue(
        default=(0, 130),
        doc="Range for the inlet flowrate of O2 [in kg/s]",
    ))

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        design_blk = self.config.design_blk

        self.initial_holdup = Var(
            within=NonNegativeReals,
            doc="Initial holdup of the tank [in tons]",
        )
        self.final_holdup = Var(
            within=NonNegativeReals,
            doc="Final holdup of the tank [in tons]",
        )
        self.lox_in = Var(
            within=NonNegativeReals,
            doc="Flowrate of LOX entering the tank [in kg/s]",
        )
        self.lox_out = Var(
            within=NonNegativeReals,
            doc="Flowrate of LOX leaving the tank [in kg/s]",
        )

        # LOX cannot be stored/withdrawn if the tank is not built
        o2_flow_range = self.config.o2_flow_range
        self.lox_in_ub_con = Constraint(
            expr=self.lox_in <= design_blk.build_tank * o2_flow_range[1]
        )
        self.lox_out_ub_con = Constraint(
            expr=self.lox_out <= design_blk.build_tank * o2_flow_range[1]
        )

        # Holdup must be greater than 10% of the tank capacity at all times
        self.holdup_constraint = Constraint(
            expr=self.final_holdup >= 0.1 * design_blk.tank_capacity,
        )

        # Mass balance over the period of one hour:
        self.mass_balance = Constraint(
            expr=self.final_holdup - self.initial_holdup == (self.lox_in - self.lox_out) * (3600 / 1000),
        )

        # Power requirement for discharging 
        self.power = Expression(
            expr=1.47054 * (self.Lox_out / 112.0865459),
            doc="Power requirement for discharging LOX [in MW]",
        )
