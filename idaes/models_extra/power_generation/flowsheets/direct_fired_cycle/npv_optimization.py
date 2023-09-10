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


from pyomo.environ import (
    ConcreteModel, 
    Var, 
    NonNegativeReals,
    Constraint,
    Expression,
    Block,
    Objective,
    maximize,
    SolverFactory,
    value,
)
from idaes.apps.grid_integration import MultiPeriodModel

from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    unit_models import (
    get_lmp_data,
    get_natural_gas_price,
    DFCDesign,
    MonoASUDesign,
    NLUDesign,
    OxygenTankDesign,
)

from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    dfc_flowsheet import (
    build_dfc_flowsheet,
    build_dfc_flowsheet_with_lox,
    build_dfc_flowsheet_with_nlu,
    append_op_costs_dfc,
    append_cashflows,
)

from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    dfc_startup_shutdown import dfc_startup_shutdown_constraints
from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    asu_startup_shutdown import asu_startup_shutdown_constraints
from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    analyze_results import _write_results
import idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    default_model_parameters as dmp


def get_linking_var_pairs(m1, m2):
    return [(m1.fs.tank.final_holdup, m2.fs.tank.initial_holdup)]


def npv_model_dfc_asu(
    dataset="MiNg_$100_CAISO_2035",
    carbon_tax=100,
    dfc_params=dmp.DFC_PARAMS,
    asu_params=dmp.MONO_ASU_PARAMS,
    cost_params=dmp.CASHFLOW_PARAMS,
):
    cost_ng = get_natural_gas_price(dataset)
    penalty = cost_params["electricity_cost"] / 1e3

    m = ConcreteModel()
    get_lmp_data(m, dataset=dataset)

    m.dfc_design = DFCDesign(model_params=dfc_params)
    m.asu_design = MonoASUDesign(model_params=asu_params)

    m.mp_model = MultiPeriodModel(
        n_time_points=365 * 24,
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
    append_cashflows(m, cost_params)

    # Declare the objective function
    m.obj = Objective(expr=m.npv, sense=maximize)

    return m


def npv_model_dfc_asu_with_lox(
    dataset="NREL",
    location="PJM-W",
    carbon_tax=150,
    dfc_params=dmp.DFC_PARAMS,
    asu_params=dmp.ASU_PARAMS,
    tank_params=dmp.O2_TANK_PARAMS,
    cost_params=dmp.CASHFLOW_PARAMS,
    solver=None,
    folder="",
):
    if dataset == "NREL":
        cost_ng = NG_PRICE_DATA[dataset][carbon_tax][location]
    else:
        cost_ng = NG_PRICE_DATA[dataset]

    penalty = cost_params["electricity_cost"] / 1e3

    m = ConcreteModel()
    get_lmp_data(m, dataset=dataset, location=location, carbon_tax=carbon_tax)

    m.dfc_design = DFCDesign(model_params=dfc_params)
    m.asu_design = MonoASUDesign(
        o2_flow_range=(10, 130), model_params=asu_params,
    )
    m.tank_design = OxygenTankDesign(
        tank_size_range=(10, 400000), model_params=tank_params,
    )

    # Currently, we are not varying the capacity of the power cycle, so fixing the value
    # of capacity. Also, we want the power cycle to be built, so fixing build_dfc. 
    m.dfc_design.build_dfc.fix(1)
    m.dfc_design.capacity.fix()

    m.mp_model = MultiPeriodModel(
        n_time_points=365 * 24,
        process_model_func=build_dfc_flowsheet_with_lox,
        linking_variable_func=get_linking_var_pairs,
        use_stochastic_build=True,
        flowsheet_options={
            "dfc_design": m.dfc_design,
            "asu_design": m.asu_design,
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
            m.mp_model.period[8760].fs.tank.holdup
        )

    # Append the overall cashflows
    append_cashflows(m, cost_params)

    # Declare the objective function
    m.obj = Objective(expr=m.npv, sense=maximize)

    # Use Gurobi solver
    if solver is None:
        solver = SolverFactory("gurobi")
        solver.options['NonConvex'] = 2
        solver.options['MIPGap'] = 0.01
        solver.options['TimeLimit'] = 7500
        solver.options['OutputFlag'] = 1

    sol = solver.solve(m, tee=True)

    _filename = folder + dataset + "_" + location + "_" + str(carbon_tax)
    _write_results(m, sol, filename=_filename, lox_withdrawal=True)

    return m, sol


def npv_model_dfc_asu_nlu(
    dataset="NREL",
    location="PJM-W",
    carbon_tax=150,
    dfc_params=dmp.DFC_PARAMS,
    asu_params=dmp.ASU_PARAMS,
    nlu_params=dmp.NLU_PARAMS,
    tank_params=dmp.O2_TANK_PARAMS,
    cost_params=dmp.CASHFLOW_PARAMS,
    solver=None,
    folder="",
):
    
    penalty = cost_params["electricity_cost"] / 1e3

    if dataset == "NREL":
        cost_ng = NG_PRICE_DATA[dataset][carbon_tax][location]
    else:
        cost_ng = NG_PRICE_DATA[dataset]

    m = ConcreteModel()
    get_lmp_data(m, dataset=dataset, location=location, carbon_tax=carbon_tax)

    m.dfc_design = DFCDesign(model_params=dfc_params)
    m.asu_design = MonoASUDesign(
        o2_flow_range=(10, 130), model_params=asu_params,
    )
    m.nlu_design = NLUDesign(
        o2_flow_range=(10, 130), model_params=nlu_params,
    )
    m.tank_design = OxygenTankDesign(
        tank_size_range=(10, 400000), model_params=tank_params,
    )

    # Currently, we are not varying the capacity of the power cycle, so fixing the value
    # of capacity. Also, we want the power cycle to be built, so fixing build_dfc. 
    m.dfc_design.build_dfc.fix(1)
    m.dfc_design.capacity.fix()

    m.mp_model = MultiPeriodModel(
        n_time_points=365 * 24,
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
            m.mp_model.period[8760].fs.tank.holdup
        )

    # Append the overall cashflows
    append_cashflows(m, cost_params)

    # Declare the objective function
    m.obj = Objective(expr=m.npv, sense=maximize)

    # Use Gurobi solver
    if solver is None:
        solver = SolverFactory("gurobi")
        solver.options['MIPGap'] = 0.01
        solver.options['TimeLimit'] = 7200
        solver.options['OutputFlag'] = 1

    sol = solver.solve(m, tee=True)

    _filename = folder + dataset + "_" + location + "_" + str(carbon_tax)
    _write_results(m, sol, filename=_filename, includes_nlu=True)

    return m, sol
