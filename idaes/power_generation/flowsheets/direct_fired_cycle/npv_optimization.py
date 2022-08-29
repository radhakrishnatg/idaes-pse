import os
import pandas as pd

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

from .unit_models import (
    get_lmp_data,
    DFCDesign,
    MonoASUDesign,
    NLUDesign,
    OxygenTankDesign,
)

from .dfc_flowsheet import (
    build_dfc_flowsheet,
    build_dfc_flowsheet_with_lox,
    build_dfc_flowsheet_with_nlu,
    append_op_costs_dfc,
    append_cashflows,
)

from .dfc_startup_shutdown import dfc_startup_shutdown_constraints
from .asu_startup_shutdown import asu_startup_shutdown_constraints


# Source: DOE ARPA - E
NG_PRICE_DATA = {
    "NREL": {
        100: {
            "CAISO": 2.26,
            "ERCOT": 2.64,
            "MISO-W": 1.69,
            "NYISO": 1.12,
            "PJM-W": 1.42
        },
        150: {
            "CAISO": 2.26,
            "ERCOT": 2.64,
            "MISO-W": 2.01,
            "NYISO": 1.14,
            "PJM-W": 1.43
        }
    },
    "Princeton": 2.94,
    "New_Princeton": 2.94
}


def _write_results(
    m,
    lox_withdrawal=False,
    includes_nlu=False,
    filename="results",
):
    # Create a directory to store alamo models
    cwd = os.getcwd()
    filename = cwd + "\\" + filename + ".xlsx"

    dfc_schedule = [m.mp_model.period[t].fs.dfc.op_mode.value for t in m.set_time]
    dfc_startup = [m.mp_model.period[t].fs.dfc.startup.value for t in m.set_time]
    dfc_shutdown = [m.mp_model.period[t].fs.dfc.shutdown.value for t in m.set_time]
    dfc_power = [m.mp_model.period[t].fs.dfc.power.value for t in m.set_time]
    asu_schedule = [m.mp_model.period[t].fs.asu.op_mode.value for t in m.set_time]
    asu_startup = [m.mp_model.period[t].fs.asu.startup.value for t in m.set_time]
    asu_shutdown = [m.mp_model.period[t].fs.asu.startup.value for t in m.set_time]
    asu_power = [value(m.mp_model.period[t].fs.asu.total_power) for t in m.set_time]
    power_to_grid = [m.mp_model.period[t].fs.power_dfc_to_grid.value for t in m.set_time]
    power_from_grid = [m.mp_model.period[t].fs.power_grid_to_asu.value for t in m.set_time]

    results = {
        "DFC_Schedule": dfc_schedule,
        "DFC_Startup": dfc_startup,
        "DFC_Shutdown": dfc_shutdown,
        "DFC_Power": dfc_power,
        "ASU_Schedule": asu_schedule,
        "ASU_Startup": asu_startup,
        "ASU_Shutdown": asu_shutdown,
        "ASU_Power": asu_power,
        "Power_to_grid": power_to_grid,
        "Power_from_grid": power_from_grid,
    }

    results_df = pd.DataFrame(results)
    results_df.to_excel(filename)


def get_linking_var_pairs(m1, m2):
    return [(m1.fs.tank.final_holdup, m2.fs.tank.initial_holdup)]


def npv_model_dfc_asu(
    dataset="NREL",
    location="PJM-W",
    carbon_tax=150,
):
    if dataset == "NREL":
        cost_ng = NG_PRICE_DATA[dataset][carbon_tax][location]
    else:
        cost_ng = NG_PRICE_DATA[dataset]

    m = ConcreteModel()
    get_lmp_data(m, dataset=dataset, location=location, carbon_tax=carbon_tax)

    m.dfc_design = DFCDesign()
    m.asu_design = MonoASUDesign()

    # Currently, we are not varying the capacity of the power cycle, so fixing the value
    # of capacity. Also, we want the power cycle to be built, so fixing build_dfc. 
    m.dfc_design.build_dfc.fix(1)
    m.dfc_design.capacity.fix()

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
        append_op_costs_dfc(m=m.mp_model.period[t], lmp=m.LMP[t], cost_ng=cost_ng)

    # Append startup and shutdown constraints for the power cycle
    m.mp_model.dfc_su_sd = Block(rule=dfc_startup_shutdown_constraints)

    # Append startup and shutdown constraints for the ASU unit
    m.mp_model.asu_su_sd = Block(rule=asu_startup_shutdown_constraints)

    # Append the overall cashflows
    append_cashflows(m)

    # Declare the objective function
    m.obj = Objective(expr=m.npv, sense=maximize)

    # Use Gurobi solver
    solver = SolverFactory("gurobi")
    solver.options['NonConvex'] = 2
    solver.options['MIPGap'] = 0.01
    solver.options['TimeLimit'] = 7500
    solver.options['OutputFlag'] = 1

    # Use BARON solver
    # solver = SolverFactory("baron")
    # solver.options["epsr"] = 0.01
    # solver.options["maxtime"] = 7500

    # Use SCIP solver
    # solver = SolverFactory("scip")

    solver.solve(m, tee=True)

    _write_results(m, filename=dataset + "_" + location + "_" + str(carbon_tax))

    return m


def npv_model_dfc_asu_nlu(
    dataset="NREL",
    location="PJM-W",
    carbon_tax=150,
):
    if dataset == "NREL":
        cost_ng = NG_PRICE_DATA[dataset][carbon_tax][location]
    else:
        cost_ng = NG_PRICE_DATA[dataset]

    m = ConcreteModel()
    get_lmp_data(m, dataset=dataset, location=location, carbon_tax=carbon_tax)

    m.dfc_design = DFCDesign()
    m.asu_design = MonoASUDesign()
    m.nlu_design = NLUDesign()
    m.tank_design = OxygenTankDesign()

    # Currently, we are not varying the capacity of the power cycle, so fixing the value
    # of capacity. Also, we want the power cycle to be built, so fixing build_dfc. 
    m.dfc_design.build_dfc.fix(1)
    m.dfc_design.capacity.fix()

    m.mp_model = MultiPeriodModel(
        n_time_points=365 * 24,
        process_model_func=build_dfc_flowsheet_with_nlu,
        linking_variable_func=None,
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
        append_op_costs_dfc(m=m.mp_model.period[t], lmp=m.LMP[t], cost_ng=cost_ng)

    # Append startup and shutdown constraints for the power cycle
    m.mp_model.dfc_su_sd = Block(rule=dfc_startup_shutdown_constraints)

    # Append startup and shutdown constraints for the ASU unit
    m.mp_model.asu_su_sd = Block(rule=asu_startup_shutdown_constraints)

    # Append the overall cashflows
    append_cashflows(m)

    # Declare the objective function
    m.obj = Objective(expr=m.npv, sense=maximize)

    # Use Gurobi solver
    solver = SolverFactory("gurobi")
    solver.options['NonConvex'] = 2
    solver.options['MIPGap'] = 0.01
    solver.options['TimeLimit'] = 7500
    solver.options['OutputFlag'] = 1

    # Use BARON solver
    # solver = SolverFactory("baron")
    # solver.options["epsr"] = 0.01
    # solver.options["maxtime"] = 7500

    # Use SCIP solver
    # solver = SolverFactory("scip")

    solver.solve(m, tee=True)

    _write_results(m, filename=dataset + "_" + location + "_" + str(carbon_tax))

    return m
