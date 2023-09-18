#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

"""
Assumptions for air separation unit:  
1. 1 kmol/hr of feed air, composition (0.7812, 0.0093, 0.2095) = (N2, Ar, O2)
2. 0.20603 kmol/hr of GOX, composition (0, 0.005, 0.995)
3. 0.63597 kmol/hr of GAN, composition (1, 0, 0)
4. 0.15 kmol/hr of WAN, composition (0.9682, 0.001852, 0.029947)
5. 0.008 kmol/hr of LAR, composition (0, 0.999, 0.001)

Molecular weight of feed: 28.964759 kg/kmol
Molecular weight of GOX : 32.03975 kg/kmol
Molecular weight of GAN : 28.02 kg/kmol
Molecular weight of WAN : 28.16128601007286 kg/kmol
Molecular weight of LAR : 39.94205 kg/kmol

This implies:
MW_feed kg air produces 0.20603 * MW_GOX kg of GOX
==> 1 kg/s of GOX requires MW_feed / (0.20603 * MW_GOX)

MW_feed kg of air produces 0.63597 * MW_GAN kg of GAN
==> MW_feed / (0.20603 * MW_GOX) kg/s of air produces 
    (0.63597 * MW_GAN) / (0.20603 * MW_GOX) of GAN

1 kg/s of Oxygen product produces:
    (0.15 * MW_WAN) / (0.20603 * MW_GOX) of WAN, and
    (0.008 * MW_LAR) / (0.20603 * MW_GOX) of LAR

Therefore, 
GAN production rate = ((0.63597 * MW_GAN) / (0.20603 * MW_GOX)) * GOX flowrate
LAR production rate = ((0.008 * MW_LAR) / (0.20603 * MW_GOX)) * GOX flowrate

"""


from importlib import resources
from pathlib import Path
import pandas as pd
import json

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
import idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    default_model_parameters as dmp

GAN_PROD_RATE = (0.63597 * 28.02) / (0.20603 * 32.03975)  # See notes above
LAR_PROD_RATE = (0.008 * 39.94205) / (0.20603 * 32.03975) # See notes above


def get_lmp_data(m, dataset="NREL"):
    """
    This function reads and appends LMP data to the model.
    """
    with resources.path(
        "idaes.models_extra.power_generation.flowsheets.direct_fired_cycle", 
        "lmp_data.csv"
    ) as p:
        path_to_file = Path(p).resolve()

    raw_data = pd.read_csv(path_to_file)
    num_days = 365

    m.set_time = RangeSet(num_days * 24)
    price_all = raw_data[dataset].tolist()

    # Set prices lower than $0.01/MWh to zero to avoid numerical issues
    price_all = [i if i > 0.01 else 0 for i in price_all]

    # LMP is divided by 1000 to formulate the objective in thousands of dollars.
    # This helps in keeping the scale of the CAPEX in O(1e6), otherwise it would be O(1e9).
    m.LMP = Param(
        m.set_time,
        initialize={t + 1: lmp / 1000 for (t, lmp) in enumerate(price_all)},
        doc="Locational Marginal Prices [in 1000$/MWh]"
    )


def get_natural_gas_price(dataset):

    with open("lmp_metadata.json") as fp:
        data = json.load(fp)

    return data[dataset]["ng_price"]


@declare_process_block_class("DFCDesign")
class DFCDesignData(SkeletonUnitModelData):
    """
    This class contains variables/parameters related to the design of the direct fired cycle
    """
    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare("model_params", ConfigValue(
        doc="Dictionary containing model parameters",
    ))

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        params = self.config.model_params
        _capacity_range = params["des_capacity_range"]  # Capacity range of DFC [MW]
        _ng_flow_coeff = params["ng_flow_coeff"]        # Coefficient relating NG flow and power [kg/s/MW]
        _o2_ng_ratio = params["o2_ng_ratio"]            # Ratio of O2 to NG flowrates [-]
        _capex = params["capex"]                        # List of coefficients of the linear surrogate
        _fom_factor = params["fom_factor"]              # Ratio of FOM and CAPEX

        self.capacity = Var(
            within=NonNegativeReals,
            bounds=(0, _capacity_range[1]),
            doc="Capacity of the power plant [in MW]",
        )

        # Define a variable that informs whether the plant needs to built or not. 
        self.build_dfc = Var(
            within=Binary,
            doc="1: Plant is built, 0: Plant is not built",
        )

        # Bound the capcity of the plant in the specified range
        self.capacity_lb_con = Constraint(expr=self.capacity >= self.build_dfc * _capacity_range[0])
        self.capacity_ub_con = Constraint(expr=self.capacity <= self.build_dfc * _capacity_range[1])

        # Compute the natural gas flowrate required at maximum capacity. 
        self.ng_flow = Expression(
            expr=_ng_flow_coeff * self.capacity,
            doc="Computes the natural flowrate required [in kg/s] at full load",
        )

        self.o2_flow = Expression(
            expr=_o2_ng_ratio * self.ng_flow,
            doc="Flowrate of oxygen [in kg/s] required at full load",
        )

        # Assumuing that the capex varies linearly with capacity
        self.capex = Expression(
            expr=_capex[0] * self.capacity + _capex[1] * self.build_dfc,
            doc="CAPEX of the power cycle [in 1000$]",
        )

        # Calculating the FOM as 3.157% of the CAPEX. The percentage value is obtained from
        # 35,641.27 (FOM) / 1,128,855 (CAPEX). The FOM is also assumed to vary linearly with the capacity.
        self.fom = Expression(
            expr=_fom_factor * self.capex,
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

        self.op_mode_capacity = Var(
            within=NonNegativeReals,
            doc="Auxiliary variable for the product of op_mode and capacity",
        )
        self.startup_capacity = Var(
            within=NonNegativeReals,
            doc="Auxiliary variable for the product of startup and capacity",
        )
        self.shutdown_capacity = Var(
            within=NonNegativeReals,
            doc="Auxiliary variable for the product of shutdown and capacity",
        )

        # Linear relaxation of the product of binary variable and capacity
        self.lin_relax_1 = Constraint(
            expr=self.op_mode_capacity + self.startup_capacity + self.shutdown_capacity <=
            design_blk.capacity
        )
        self.lin_relax_2 = Constraint(
            expr=self.op_mode_capacity + self.startup_capacity + self.shutdown_capacity >=
            design_blk.capacity - design_blk.capacity.ub * (
                design_blk.build_dfc - self.op_mode - self.startup - self.shutdown
            )
        )
        self.lin_relax_3 = Constraint(
            expr=self.op_mode_capacity <= self.op_mode * design_blk.capacity.ub
        )
        self.lin_relax_4 = Constraint(
            expr=self.startup_capacity <= self.startup * design_blk.capacity.ub
        )
        self.lin_relax_5 = Constraint(
            expr=self.shutdown_capacity <= self.shutdown * design_blk.capacity.ub
        )

        # Power production as a function of natural gas flowrate
        # We construct a surrogate model of the form (1 - norm_power) = m * (1 - norm_ng_flow).
        # This way, the natural gas requirement is exact at full load.
        # Rearranging the equation yields norm_power = m * norm_ng_flow + 1-m 
        # Next, we replace norm_power = power / capacity and norm_ng_flow = ng_flow / max_ng_flow,
        # substitute max_ng_flow in terms of capacity and rearrange the equation.

        params = design_blk.config.model_params
        _operating_range = params["op_capacity_range"]
        _ng_flow_coeff = params["ng_flow_coeff"]
        _perf_curve = params["op_curve_coeff"]
        _o2_ng_ratio = params["o2_ng_ratio"]
        _const_vom_coeff = params["const_vom_coeff"]
        _var_vom_coeff = params["var_vom_coeff"]
        _co2_emission = params["co2_emission_rate"]
        _co2_capture = params["co2_capture_rate"]

        self.power_production = Constraint(
            expr=self.power == (_perf_curve[0] / _ng_flow_coeff) * self.ng_flow + 
            _perf_curve[1] * self.op_mode_capacity,
        )

        # Ensure that power production is within P_min and P_max
        self.power_production_lb = Constraint(
            expr=_operating_range[0] * self.op_mode_capacity <= self.power,
        )

        self.power_production_ub = Constraint(
            expr=self.power <= _operating_range[1] * self.op_mode_capacity,
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
            expr=_o2_ng_ratio * self.total_ng_flow,
            doc="Flowrate of oxygen required [in kg/s]",
        )

        self.non_fuel_vom = Expression(
            expr=(_var_vom_coeff * self.power + 
                  _const_vom_coeff[0] * self.op_mode_capacity + 
                  _const_vom_coeff[1] * self.op_mode),
            doc="Non-fuel VOM cost [in $1000/hr]",
        )

        self.co2_emissions = Expression(
            expr=_co2_emission * self.total_ng_flow * (1 - _co2_capture),
            doc="Net CO2 emissions [in kg/s]",
        )

    def change_shutdown_to_continuous(self):
        self.shutdown.domain = NonNegativeReals
        self.shutdown.setub(1)

    def change_shutdown_to_binary(self):
        self.shutdown.domain = Binary


@declare_process_block_class("MonoASUDesign")
class MonoASUDesignData(SkeletonUnitModelData):
    """
    This class builds the model for the design of the ASU. Here, we treat the
    system of ASUs as a monolithic ASU.
    """

    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare("model_params", ConfigValue(
        doc="Dictionary containing model parameters",
    ))

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        params = self.config.model_params
        _o2_flow_range = params["des_capacity_range"]
        _power_req = params["power_req_coeff"]
        _capex = params["capex"]
        _fom_factor = params["fom_factor"]

        self.max_o2_flow = Var(
            within=NonNegativeReals,
            bounds=(0, _o2_flow_range[1]),
            doc="Maximum flowrate of O2 the ASU can produce [in kg/s]",
        )
        self.build_asu = Var(
            within=Binary,
            doc="1: ASU is built, 0: ASU is not built",
        )

        # Bound the capcity of the plant in the specified range
        self.o2_flow_lb_con = Constraint(expr=self.max_o2_flow >= self.build_asu * _o2_flow_range[0])
        self.o2_flow_ub_con = Constraint(expr=self.max_o2_flow <= self.build_asu * _o2_flow_range[1])

        # Relation between the flowrate and the power requirement
        self.max_power = Expression(
            expr=_power_req * self.max_o2_flow,
            doc="Power requirement at maximum capacity [in MW]",
        )

        # Assuming that the capex of the ASU varies linearly with size
        self.capex = Expression(
            expr=_capex[0] * self.max_o2_flow + _capex[1] * self.build_asu,
            doc="CAPEX of the ASU unit [in 1000$]",
        )

        # Calculating the FOM as 3.157% of the CAPEX. The percentage value is obtained from
        # 17223.73 (FOM) / 545,522 (CAPEX). The FOM is also assumed to vary linearly with the capacity. 
        self.fom = Expression(
            expr=_fom_factor * self.capex,
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

        self.op_mode_o2_flow = Var(
            within=NonNegativeReals,
            doc="Auxiliary variable for linearizing the product of op_mode and max_o2_flow"
        )
        self.startup_o2_flow = Var(
            within=NonNegativeReals,
            doc="Auxiliary variable for linearizing the product of startup and max_o2_flow"
        )
        self.shutdown_o2_flow = Var(
            within=NonNegativeReals,
            doc="Auxiliary variable for linearizing the product of shutdown and max_o2_flow"
        )

        # Linear relaxation of op_mode * max_o2_flow
        self.lin_relax_1 = Constraint(
            expr=self.op_mode_o2_flow + self.startup_o2_flow + self.shutdown_o2_flow <=
            design_blk.max_o2_flow
        )
        self.lin_relax_2 = Constraint(
            expr=self.op_mode_o2_flow + self.startup_o2_flow + self.shutdown_o2_flow >= 
            design_blk.max_o2_flow - design_blk.max_o2_flow.ub * (
                design_blk.build_asu - self.op_mode - self.startup - self.shutdown
            )
        )
        self.lin_relax_3 = Constraint(
            expr=self.op_mode_o2_flow <= self.op_mode * design_blk.max_o2_flow.ub
        )
        self.lin_relax_4 = Constraint(
            expr=self.startup_o2_flow <= self.startup * design_blk.max_o2_flow.ub
        )
        self.lin_relax_5 = Constraint(
            expr=self.shutdown_o2_flow <= self.shutdown * design_blk.max_o2_flow.ub
        )

        # Power production as a function of natural gas flowrate
        # We construct a surrogate model of the form (1 - norm_power) = m * (1 - norm_o2_flow).
        # This way, the power requirement is exact at full capacity.
        # Rearranging the equation yields norm_power = m * norm_o2_flow + 1-m 

        params = design_blk.config.model_params
        _operating_range = params["op_capacity_range"]
        _power_req = params["power_req_coeff"]
        _perf_curve = params["op_curve_coeff"]
        _var_vom = params["var_vom_coeff"]
        _const_vom = params["const_vom_coeff"]
        _nitrogen_price = params["nitrogen_price"]
        _argon_price = params["argon_price"]

        self.power_requirement = Constraint(
            expr=(1 / _power_req) * self.power == 
            _perf_curve[0] * self.o2_flow + _perf_curve[1] * self.op_mode_o2_flow, 
        )

        # Ensure that the normalized oxygen flowrate is within the admissible operating range
        self.o2_flow_lb = Constraint(
            expr=_operating_range[0] * self.op_mode_o2_flow <= self.o2_flow,
        )

        self.o2_flow_ub = Constraint(
            expr=self.o2_flow <= _operating_range[1] * self.op_mode_o2_flow,
        )

        self.su_sd_power = Var(
            within=NonNegativeReals,
            doc="Power requirement during startup and shutdown [in MW]",
        )

        self.total_power = Expression(
            expr=self.power + self.su_sd_power,
            doc="Total power required by the ASU [in MW]",
        )

        self.non_fuel_vom = Expression(
            expr=(_var_vom * self.power + _const_vom[0] * self.op_mode_o2_flow 
                  + _const_vom[1] * self.op_mode),
            doc="Non-electricity VOM cost [in $1000/hr]",
        )

        self.nitrogen_revenue = Expression(
            expr=GAN_PROD_RATE * self.o2_flow * 3600 * _nitrogen_price,
            doc="Revenue from nitrogen market [$/hr]",
        )

        self.argon_revenue = Expression(
            expr=LAR_PROD_RATE * self.o2_flow * 3600 * _argon_price,
            doc="Revenue from argon market [$/hr]",
        )

    def change_shutdown_to_continuous(self):
        self.shutdown.domain = NonNegativeReals
        self.shutdown.setub(1)

    def change_shutdown_to_binary(self):
        self.shutdown.domain = Binary


@declare_process_block_class("NLUDesign")
class NLUDesignData(SkeletonUnitModelData):
    """
    This class builds the model for the design of the NLU.
    """

    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare("model_params", ConfigValue(
        doc="Dictionary containing NLU model parameters",
    ))

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        params = self.config.model_params
        _o2_flow_range = params["des_capacity_range"]
        _power_req = params["power_req_coeff"]
        _capex = params["capex"]
        _fom_factor = params["fom_factor"]

        self.max_o2_flow = Var(
            within=NonNegativeReals,
            bounds=(0, _o2_flow_range[1]),
            doc="Maximum flowrate of O2 the ASU can produce [in kg/s]",
        )
        self.build_nlu = Var(
            within=Binary,
            doc="1: NLU is built, 0: NLU is not built",
        )

        # Bound the capcity of the plant in the specified range
        self.o2_flow_lb_con = Constraint(expr=self.max_o2_flow >= self.build_nlu * _o2_flow_range[0])
        self.o2_flow_ub_con = Constraint(expr=self.max_o2_flow <= self.build_nlu * _o2_flow_range[1])

        # Relation between the flowrate and the power requirement
        self.max_power = Expression(
            expr=_power_req * self.max_o2_flow,
            doc="Power requirement at maximum capacity [in MW]",
        )

        # Assuming that the capex of the NLU varies linearly with size
        self.capex = Expression(
            expr=_capex[0] * self.max_o2_flow + _capex[1] * self.build_nlu,
            doc="CAPEX of the ASU unit [in 1000$]",
        )

        # Assuming that the FOM of NLU is 3.157% of CAPEX. This number is calculated as
        # 41,715 (FOM) / 171,189 (CAPEX) 
        self.fom = Expression(
            expr=_fom_factor * self.capex,
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

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        design_blk = self.config.design_blk
        params = design_blk.config.model_params
        _power_req = params["power_req_coeff"]
        _operating_range = params["op_capacity_range"]
        _var_vom = params["var_vom_coeff"]
        _const_vom = params["const_vom_coeff"]

        # Declare variables
        self.power = Var(
            within=NonNegativeReals,
            doc="Power required to liquify oxygen [in MW]",
        )
        self.o2_flow = Var(
            within=NonNegativeReals,
            doc="Flowrate of liquified oxygen [in kg/s]",
        )

        # Note: Assuming that the NLU does not have any startup/shutdown constraints.
        # If it does, we need to define startup and shutdown binary variables like before.
        self.op_mode = Var(
            within=Binary,
            doc="1: In Operation, 0: Shutdown",
        )
        self.op_mode_o2_flow = Var(
            within=NonNegativeReals,
            doc="Auxiliary variable for linearizing op_mode and max_o2_flow",
        )
        self.mccor_conv = Constraint(
            expr=self.op_mode_o2_flow >= design_blk.max_o2_flow 
            + self.op_mode * design_blk.max_o2_flow.ub 
            - design_blk.max_o2_flow.ub * design_blk.build_nlu,
        )
        self.mccor_conc_1 = Constraint(
            expr=self.op_mode_o2_flow <= design_blk.max_o2_flow,
        )
        self.mccor_conc_2 = Constraint(
            expr=self.op_mode_o2_flow <= self.op_mode * design_blk.max_o2_flow.ub
        )

        # Power production as a function of natural gas flowrate
        self.power_requirement = Constraint(
            expr=self.power == _power_req * self.o2_flow,
        )

        # Ensure that the normalized oxygen flowrate is within the operating_range
        self.o2_flow_lb = Constraint(
            expr=_operating_range[0] * self.op_mode_o2_flow <= self.o2_flow,
        )

        self.o2_flow_ub = Constraint(
            expr=self.o2_flow <= _operating_range[1] * self.op_mode_o2_flow,
        )

        self.non_fuel_vom = Expression(
            expr=(_var_vom * self.power + _const_vom[0] * self.op_mode_o2_flow 
                  + _const_vom[1] * self.op_mode),
            doc="Non-electricity VOM cost [in $1000/hr]",
        )


@declare_process_block_class("OxygenTankDesign")
class OxygenTankDesignData(SkeletonUnitModelData):
    """
    This class builds the model for the design of the tank.
    """

    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare("model_params", ConfigValue(
        doc="Dictionary containing model parameters"
    ))

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        params = self.config.model_params
        _tank_size_range = params["des_capacity_range"]
        _fom_factor = params["fom_factor"]
        _capex = params["capex"]

        self.tank_capacity = Var(
            within=NonNegativeReals,
            bounds=(0, _tank_size_range[1]),
            doc="Maximum amount of oxygen the tank can hold [in tons]",
        )

        self.build_tank = Var(
            within=Binary,
            doc="1: The storage tank is built, 0: Otherwise",
        )

        # Bound the capcity of the tank in the specified range
        self.tank_capacity_lb_con = Constraint(
            expr=self.tank_capacity >= self.build_tank * _tank_size_range[0]
        )
        self.tank_capacity_ub_con = Constraint(
            expr=self.tank_capacity <= self.build_tank * _tank_size_range[1]
        )

        self.capex = Expression(
            expr=_capex[0] * self.tank_capacity + _capex[1] * self.build_tank,
            doc="CAPEX of the oxygen storage tank [in 1000$]",
        )

        self.fom = Expression(
            expr=_fom_factor * self.capex,
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

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        design_blk = self.config.design_blk
        params = design_blk.config.model_params
        _min_holdup = params["min_holdup"]
        _power_req = params["power_req_coeff"]
        _max_o2_flow = 130  # Maximum flowrate of Oxygen entering the tank 

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
        self.lox_in_ub_con = Constraint(
            expr=self.lox_in <= design_blk.build_tank * _max_o2_flow
        )
        self.lox_out_ub_con = Constraint(
            expr=self.lox_out <= design_blk.build_tank * _max_o2_flow
        )

        # Holdup must be greater than 10% of the tank capacity at all times
        self.holdup_constraint = Constraint(
            expr=self.final_holdup >= _min_holdup * design_blk.tank_capacity,
        )

        # Holdup cannot exceed the capacity of the tank
        self.holdup_ub = Constraint(
            expr=self.final_holdup <= design_blk.tank_capacity
        )

        # Mass balance over the period of one hour:
        self.mass_balance = Constraint(
            expr=self.final_holdup - self.initial_holdup == (self.lox_in - self.lox_out) * (3600 / 1000),
        )

        # Power requirement for discharging 
        self.power = Expression(
            expr=_power_req * self.lox_out,
            doc="Power requirement for discharging LOX [in MW]",
        )
