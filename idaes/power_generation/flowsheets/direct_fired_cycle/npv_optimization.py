import os
import pandas as pd
import json

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
    sol,
    lox_withdrawal=False,
    includes_nlu=False,
    filename="results",
):
    # Create a directory to store NPV optimization results
    cwd = os.getcwd()
    _filename = cwd + "\\npv_results\\" + filename + ".xlsx"

    if not os.path.exists(cwd + "\\npv_results"):
        os.mkdir(os.path.join(cwd, "npv_results"))

    set_flowsheets = m.mp_model.period[:].fs

    results = {
        "LMP [$/MWh]": [m.LMP[t] * 1000 for t in m.set_time],
        "DFC_Schedule": [value(fs.dfc.op_mode) for fs in set_flowsheets],
        "DFC_Startup": [value(fs.dfc.startup) for fs in set_flowsheets],
        "DFC_Shutdown": [value(fs.dfc.shutdown) for fs in set_flowsheets],
        "DFC_Power": [value(fs.dfc.power) for fs in set_flowsheets],
        "DFC_Ng_Flow": [value(fs.dfc.total_ng_flow) for fs in set_flowsheets],
        "DFC_O2_Flow": [value(fs.dfc.o2_flow) for fs in set_flowsheets],

        "ASU_Schedule": [value(fs.asu.op_mode) for fs in set_flowsheets],
        "ASU_Startup": [value(fs.asu.startup) for fs in set_flowsheets],
        "ASU_Shutdown": [value(fs.asu.shutdown) for fs in set_flowsheets],
        "ASU_Power": [value(fs.asu.total_power) for fs in set_flowsheets],
        "ASU_O2_Flow": [value(fs.asu.o2_flow) for fs in set_flowsheets],

        "Power_to_grid": [value(fs.power_dfc_to_grid) for fs in set_flowsheets],
        "Power_grid_to_asu": [value(fs.power_grid_to_asu) for fs in set_flowsheets],
        "Power_dfc_to_asu": [value(fs.power_dfc_to_asu) for fs in set_flowsheets],

        "Oxygen_asu_to_dfc": [value(fs.oxygen_asu_to_dfc) for fs in set_flowsheets],
        "Oxygen_asu_to_vent": [value(fs.oxygen_asu_to_vent) for fs in set_flowsheets],
    }

    if includes_nlu:
        results["NLU_Schedule"] = [value(fs.nlu.op_mode) for fs in set_flowsheets]
        results["NLU_O2_Flow"] = [value(fs.nlu.o2_flow) for fs in set_flowsheets]
        results["NLU_Power"] = [value(fs.nlu.power) for fs in set_flowsheets]

        results["Tank_init_holdup"] = [value(fs.tank.initial_holdup) for fs in set_flowsheets]
        results["Tank_final_holdup"] = [value(fs.tank.final_holdup) for fs in set_flowsheets]
        results["Tank_lox_in"] = [value(fs.tank.lox_in) for fs in set_flowsheets]
        results["Tank_lox_out"] = [value(fs.tank.lox_out) for fs in set_flowsheets]

        results["Power_dfc_to_nlu"] = [value(fs.power_dfc_to_nlu) for fs in set_flowsheets]
        results["Power_dfc_to_tank"] = [value(fs.power_dfc_to_Tank) for fs in set_flowsheets]
        results["Power_grid_to_nlu"] = [value(fs.power_grid_to_nlu) for fs in set_flowsheets]
        results["Power_grid_to_tank"] = [value(fs.power_grid_to_tank) for fs in set_flowsheets]

    if lox_withdrawal:
        results["GOx_fraction"] = [value(fs.gox_fraction) for fs in set_flowsheets]

        results["Tank_init_holdup"] = [value(fs.tank.initial_holdup) for fs in set_flowsheets]
        results["Tank_final_holdup"] = [value(fs.tank.final_holdup) for fs in set_flowsheets]
        results["Tank_lox_in"] = [value(fs.tank.lox_in) for fs in set_flowsheets]
        results["Tank_lox_out"] = [value(fs.tank.lox_out) for fs in set_flowsheets]

        results["Power_dfc_to_tank"] = [value(fs.power_dfc_to_Tank) for fs in set_flowsheets]
        results["Power_grid_to_tank"] = [value(fs.power_grid_to_tank) for fs in set_flowsheets]

    results_df = pd.DataFrame(results)
    results_df.to_excel(_filename)

    lower_bnd = sol["Problem"][0]["Lower bound"]
    upper_bnd = sol["Problem"][0]["Upper bound"]

    solution = {
        "Lower bound": lower_bnd,
        "Upper bound": upper_bnd,
        "Gap": ((upper_bnd - lower_bnd) / lower_bnd) * 100,
        "Wall time": sol["Solver"][0]["Wall time"],
        "Status": sol["Solver"][0]["Status"],
        "Termination message": sol["Solver"][0]["Termination message"],

        "NPV": value(m.obj) / 1000,
        "DFC_Capacity": m.dfc_design.capacity.value,
        "ASU_Capacity": m.asu_design.max_o2_flow.value,
        "NLU_Capacity": (m.nlu_design.max_o2_flow.value if hasattr(m, "nlu_design") else "N/A"),
        "Tank_Capacity": (m.tank_design.tank_capacity.value 
                          if hasattr(m, "tank_design") else "N/A"),
        "DFC Startups": sum(results["DFC_Startup"]),
        "DFC Shutdowns": sum(results["DFC_Shutdown"]),
        "ASU Startups": sum(results["ASU_Startup"]),
        "ASU Shutdowns": sum(results["ASU_Shutdown"]),
        "Total Power Produced": sum(results["DFC_Power"]),
        "Total Power Sold": sum(results["Power_to_grid"]),
        "Total CAPEX": m.CAPEX.value / 1000,
        "Total FOM": m.FOM.value / 1000,
        "Total Net Profit": m.NET_PROFIT.value / 1000,

        "Num Vars": sol["Problem"][0]["Number of variables"],
        "Num Bin Vars": sol["Problem"][0]["Number of binary variables"],
        "Num constraints": sol["Problem"][0]["Number of constraints"],
    }

    with open(cwd + "\\npv_results\\" + filename + ".json", "w") as fp:
        json.dump(solution, fp, indent=4)


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
        append_op_costs_dfc(m=m.mp_model.period[t], lmp=m.LMP[t], cost_ng=cost_ng, carbon_price=carbon_tax / 1000)

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
    # solver.options["limits/gap"] = 0.01
    # solver.options["limits/time"] = 7500
    # solver.options["lp/threads"] = 16

    sol = solver.solve(m, tee=True)

    _write_results(m, sol, filename=dataset + "_" + location + "_" + str(carbon_tax))

    return m, sol


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
        append_op_costs_dfc(m=m.mp_model.period[t], lmp=m.LMP[t], cost_ng=cost_ng, carbon_price=carbon_tax / 1000)

    # Append startup and shutdown constraints for the power cycle
    m.mp_model.dfc_su_sd = Block(rule=dfc_startup_shutdown_constraints)

    # Append startup and shutdown constraints for the ASU unit
    m.mp_model.asu_su_sd = Block(rule=asu_startup_shutdown_constraints)

    # Set the initial holdup of the tank
    m.mp_model.initial_tank_level = Constraint(
        expr=m.mp_model.period[1].fs.tank.initial_holdup == 0.1 * m.tank_design.tank_capacity
    )

    # Append the overall cashflows
    append_cashflows(m)

    # Declare the objective function
    m.obj = Objective(expr=m.npv, sense=maximize)

    # Use Gurobi solver
    solver = SolverFactory("gurobi")
    solver.options['NonConvex'] = 2
    solver.options['MIPGap'] = 0.01
    solver.options['TimeLimit'] = 15000
    solver.options['OutputFlag'] = 1

    # Use BARON solver
    # solver = SolverFactory("baron")
    # solver.options["epsr"] = 0.01
    # solver.options["maxtime"] = 7500

    # Use SCIP solver
    # solver = SolverFactory("scip")
    # solver.options["limits/gap"] = 0.01
    # solver.options["limits/time"] = 7500
    # solver.options["lp/threads"] = 16

    sol = solver.solve(m, tee=True)

    _write_results(m, sol, filename=dataset + "_" + location + "_" + str(carbon_tax), includes_nlu=True)

    return m, sol
