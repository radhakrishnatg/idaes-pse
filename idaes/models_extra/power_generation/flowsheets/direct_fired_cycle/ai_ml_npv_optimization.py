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

from importlib import resources
from pathlib import Path
import pandas as pd

from pyomo.environ import (
    ConcreteModel, 
    Var,
    Expression,
    NonNegativeReals,
    RangeSet,
    Param,
    Constraint,
    Block,
    Objective,
    maximize,
)
from idaes.apps.grid_integration import MultiPeriodModel

from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    unit_models import (
    get_natural_gas_price,
    DFCDesign,
    MonoASUDesign,
    NLUDesign,
    OxygenTankDesign,
)

from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    dfc_flowsheet import (
    build_dfc_flowsheet,
    build_dfc_flowsheet_with_nlu,
    append_op_costs_dfc,
)

from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    dfc_startup_shutdown import dfc_startup_shutdown_constraints
from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    asu_startup_shutdown import asu_startup_shutdown_constraints
from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    analyze_results import _write_results
import idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    default_model_parameters as dmp


def get_lmp_data(m, num_time_steps=7*24, dataset="NREL"):
    """
    This function reads and appends LMP data to the model.
    """
    with resources.path(
        "idaes.models_extra.power_generation.flowsheets.direct_fired_cycle", 
        "lmp_data.csv"
    ) as p:
        path_to_file = Path(p).resolve()

    raw_data = pd.read_csv(path_to_file)

    m.set_time = RangeSet(num_time_steps)
    price_all = raw_data[dataset].tolist()

    # Set prices lower than $0.01/MWh to zero to avoid numerical issues
    price_all = [i if i > 0.01 else 0 for i in price_all]

    # LMP is divided by 1000 to formulate the objective in thousands of dollars.
    # This helps in keeping the scale of the CAPEX in O(1e6), otherwise it would be O(1e9).
    m.LMP = Param(
        m.set_time,
        initialize={t: price_all[t - 1] / 1000 for t in m.set_time},
        doc="Locational Marginal Prices [in 1000$/MWh]"
    )


def append_cashflows(m, num_time_steps, params):
    plant_life = params["plant_life"]
    tax_rate = params["tax_rate"]
    fcr = params["FCR"]

    m.CAPEX = Var(
        within=NonNegativeReals,
        doc="Total CAPEX of the plant [in $1000]",
    )
    m.FOM = Var(
        within=NonNegativeReals,
        doc="Total fixed O&M costs per year [in $1000]",
    )
    m.DEPRECIATION = Var(
        within=NonNegativeReals,
        doc="Depreciation value per year [in $1000]"
    )
    m.CORP_TAX = Var(
        within=NonNegativeReals,
        doc="Net corporate tax per year [in $1000]",
    )
    m.NET_PROFIT = Var(
        doc="Net profit per year [in $1000]",
    )

    m.capex_calculation = Constraint(
        expr=m.CAPEX == m.dfc_design.capex + m.asu_design.capex + 
        (m.nlu_design.capex if hasattr(m, "nlu_design") else 0) +
        (m.tank_design.capex if hasattr(m, "tank_design") else 0),
        doc="Calculates the total CAPEX [in $1000]",
    )
    m.fom_calculation = Constraint(
        expr=m.FOM == m.dfc_design.fom + m.asu_design.fom +
        (m.nlu_design.fom if hasattr(m, "nlu_design") else 0) +
        (m.tank_design.fom if hasattr(m, "tank_design") else 0),
        doc="Calculates the total FOM [in $1000]",
    )
    m.depreciation_calculation = Constraint(
        expr=m.DEPRECIATION == m.CAPEX / plant_life,
        doc="Straight line depreciation with zero salvage value [in $1000]",
    )
    
    # Tax = max{0, Formula below}. We will relax it to Tax >= max{0, Formula below}
    # The inequality will be binding for the optimal solution. 
    m.corp_tax_calculation = Constraint(
        expr=m.CORP_TAX >= tax_rate * (
            (8760 / num_time_steps) *  sum(m.mp_model.period[t].fs.net_cash_flow for t in m.set_time)
            - m.FOM - m.DEPRECIATION
        ),
        doc="Calculates the total corporate tax [in $1000]",
    )

    m.net_profit_calculation = Constraint(
        expr=m.NET_PROFIT == 
        + (8760 / num_time_steps) * sum(m.mp_model.period[t].fs.net_cash_flow for t in m.set_time)
        - m.FOM - m.CORP_TAX,
        doc="Calculate the net profit generated in a year [in $1000]",
    )

    # NPV Calculation
    m.npv = Expression(
        expr=m.NET_PROFIT - fcr * m.CAPEX,
        doc="Calculates the net present value [in $1000]",
    )


def get_linking_var_pairs(m1, m2):
    return [(m1.fs.tank.final_holdup, m2.fs.tank.initial_holdup)]


def npv_model_dfc_asu(
    num_time_steps=7 * 24,
    dataset="MiNg_$100_CAISO_2035",
    carbon_tax=100,
    dfc_params=dmp.DFC_PARAMS,
    asu_params=dmp.MONO_ASU_PARAMS,
    cost_params=dmp.CASHFLOW_PARAMS,
):
    cost_ng = get_natural_gas_price(dataset)
    penalty = cost_params["electricity_cost"] / 1e3

    m = ConcreteModel()
    get_lmp_data(m, num_time_steps=num_time_steps, dataset=dataset)

    m.dfc_design = DFCDesign(model_params=dfc_params)
    m.asu_design = MonoASUDesign(model_params=asu_params)

    m.mp_model = MultiPeriodModel(
        n_time_points=num_time_steps,
        process_model_func=build_dfc_flowsheet,
        linking_variable_func=None,
        use_stochastic_build=True,
        flowsheet_options={
            "dfc_design": m.dfc_design,
            "asu_design": m.asu_design,
        },
    )

    # Append cashflows at each hour
    for t in m.set_time:
        append_op_costs_dfc(
            m=m.mp_model.period[t], 
            lmp=m.LMP[t], 
            penalty=penalty, 
            cost_ng=cost_ng, 
            carbon_price=carbon_tax / 1e3,
        )

    # Append startup and shutdown constraints for the power cycle
    m.mp_model.dfc_su_sd = Block(rule=dfc_startup_shutdown_constraints)

    # Append startup and shutdown constraints for the ASU unit
    m.mp_model.asu_su_sd = Block(rule=asu_startup_shutdown_constraints)

    # Append the overall cashflows
    append_cashflows(m, num_time_steps, cost_params)

    # Declare the objective function
    m.obj = Objective(expr=m.npv, sense=maximize)

    return m


def npv_model_dfc_asu_nlu(
    num_time_steps=7 * 24,
    dataset="MiNg_$100_CAISO_2035",
    carbon_tax=100,
    dfc_params=dmp.DFC_PARAMS,
    asu_params=dmp.MONO_ASU_PARAMS,
    nlu_params=dmp.NLU_PARAMS,
    tank_params=dmp.O2_TANK_PARAMS,
    cost_params=dmp.CASHFLOW_PARAMS,
):
    
    cost_ng = get_natural_gas_price(dataset)
    penalty = cost_params["electricity_cost"] / 1e3

    m = ConcreteModel()
    get_lmp_data(m, num_time_steps=num_time_steps, dataset=dataset)

    m.dfc_design = DFCDesign(model_params=dfc_params)
    m.asu_design = MonoASUDesign(model_params=asu_params)
    m.nlu_design = NLUDesign(model_params=nlu_params)
    m.tank_design = OxygenTankDesign(model_params=tank_params)

    m.mp_model = MultiPeriodModel(
        n_time_points=num_time_steps,
        process_model_func=build_dfc_flowsheet_with_nlu,
        linking_variable_func=get_linking_var_pairs,
        use_stochastic_build=True,
        flowsheet_options={
            "dfc_design": m.dfc_design,
            "asu_design": m.asu_design,
            "nlu_design": m.nlu_design,
            "tank_design": m.tank_design,
        },
    )

    # Append cashflows at each hour
    for t in m.set_time:
        append_op_costs_dfc(
            m=m.mp_model.period[t], 
            lmp=m.LMP[t], 
            penalty=penalty, 
            cost_ng=cost_ng, 
            carbon_price=carbon_tax / 1e3,
        )

    # Append startup and shutdown constraints for the power cycle
    m.mp_model.dfc_su_sd = Block(rule=dfc_startup_shutdown_constraints)

    # Append startup and shutdown constraints for the ASU unit
    m.mp_model.asu_su_sd = Block(rule=asu_startup_shutdown_constraints)

    # Set the initial holdup of the tank
    if tank_params["tank_constraint"] == "initial_holdup":
        m.mp_model.initial_tank_level = Constraint(
            expr=m.mp_model.period[1].fs.tank.initial_holdup == 
            tank_params["min_holdup"] * m.tank_design.tank_capacity
        )

    elif tank_params["tank_constraint"] == "periodic":
        m.mp_model.periodic_constraint = Constraint(
            expr=m.mp_model.period[1].fs.tank.initial_holdup == 
            m.mp_model.period[num_time_steps].fs.tank.final_holdup
        )

    # Append the overall cashflows
    append_cashflows(m, num_time_steps, cost_params)

    # Declare the objective function
    m.obj = Objective(expr=m.npv, sense=maximize)

    return m
