import os
import copy

from pyomo.environ import SolverFactory, Constraint
from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle.npv_optimization import (
    npv_model_dfc_asu,
    npv_model_dfc_asu_with_lox,
    npv_model_dfc_asu_nlu,
)
import idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    default_model_parameters as dmp
from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle. \
    analyze_results import _write_results

DFC_PARAMS = dmp.DFC_PARAMS
ASU_PARAMS = dmp.MONO_ASU_PARAMS
TANK_PARAMS = dmp.O2_TANK_PARAMS
NLU_PARAMS = dmp.NLU_PARAMS
CASHFLOW_PARAMS = dmp.CASHFLOW_PARAMS


price_signal_list = {
    1: {"dataset": "MiNg_$100_CAISO_2035", "carbon_tax": 100},
    2: {"dataset": "MiNg_$100_ERCOT_2035", "carbon_tax": 100},
    3: {"dataset": "MiNg_$100_MISO-W_2035", "carbon_tax": 100},
    4: {"dataset": "MiNg_$100_NYISO_2035", "carbon_tax": 100},
    5: {"dataset": "MiNg_$100_PJM-W_2035", "carbon_tax": 100},
    6: {"dataset": "MiNg_$150_CAISO_2035", "carbon_tax": 150},
    7: {"dataset": "MiNg_$150_ERCOT_2035", "carbon_tax": 150},
    8: {"dataset": "MiNg_$150_MISO-W_2035", "carbon_tax": 150},
    9: {"dataset": "MiNg_$150_NYISO_2035", "carbon_tax": 150},
    10: {"dataset": "MiNg_$150_PJM-W_2035", "carbon_tax": 150},

    "V1": {"dataset": "EXEMPLAR_LARGE", "carbon_tax": 0},
    "V2": {"dataset": "EXEMPLAR_MEDIUM", "carbon_tax": 0},
    "V3": {"dataset": "EXEMPLAR_SMALL", "carbon_tax": 0},
}

solver = SolverFactory("gurobi")
solver.options['MIPGap'] = 0.01
solver.options['TimeLimit'] = 7500
solver.options['OutputFlag'] = 1


def single_case(dataset="MiNg_$100_CAISO_2035", carbon_tax=100):
    dfc_params = copy.deepcopy(DFC_PARAMS)
    asu_params = copy.deepcopy(ASU_PARAMS)
    cost_params = copy.deepcopy(CASHFLOW_PARAMS)

    # Modify any default parameters if needed
    # dfc_params["capex"] /= 2

    cwd = os.getcwd()
    if not os.path.exists(cwd + "\\single_case"):
        os.mkdir(os.path.join(cwd, "single_case"))

    m = npv_model_dfc_asu(
        dataset=dataset,
        carbon_tax=carbon_tax,
        dfc_params=dfc_params,
        asu_params=asu_params,
        cost_params=cost_params,
    )

    # Enforce that the DFC is built and fix its capacity
    m.dfc_design.build_dfc.fix(1)
    m.dfc_design.capacity.fix(777.685)

    # Solve the model
    sol = solver.solve(m, tee=True)

    # Write results to file
    _filename = "single_case/" + dataset
    _write_results(m, sol, filename=_filename) 

    return m


def validation_runs():
    dfc_params = copy.deepcopy(DFC_PARAMS)
    asu_params = copy.deepcopy(ASU_PARAMS)
    cost_params = copy.deepcopy(CASHFLOW_PARAMS)

    # Note: Capex costs in the default model parameters only include
    # plant costs. It does not include owner's costs.
    # Owners costs are roughly 21.3496% of the Total plant cost, so we multiply
    # TPC with 1.213496 to obtain the Total Owner's cost
    # Then, we multiply it with 1.093 to obtain the Total as-spent cost
    # Annualization factor (FCR) is 0.0707
    # Depreciation is included in FRC, so we set corporate tax rate to zero

    # Modify CAPEX coefficients here
    cost_params["FCR"] = 0.0707 * 1.093 * 1.213496
    cost_params["tax_rate"] = 0

    cwd = os.getcwd()
    if not os.path.exists(cwd + "\\validation"):
        os.mkdir(os.path.join(cwd, "validation"))

    dfc_size = {"V1": 777.685, "V2": 422.049, "V3": 176.032}

    for ps in ["V1", "V2", "V3"]:
        m = npv_model_dfc_asu(
            dataset=price_signal_list[ps]["dataset"],
            carbon_tax=price_signal_list[ps]["carbon_tax"],
            dfc_params=dfc_params,
            asu_params=asu_params,
            cost_params=cost_params,
        )

        # Add capacity factor constraint
        m.cap_fac_cons = Constraint(
            expr=sum(m.mp_model.period[t].fs.dfc.op_mode for t in m.mp_model.set_time) == 7446
        )

        m.cap_fac_cons_1 = Constraint(
            expr=sum(m.mp_model.period[t].fs.dfc.op_mode for t in range(1, 7447)) == 7446
        )

        # Enforce that the DFC is built and fix its capacity
        m.dfc_design.build_dfc.fix(1)
        m.dfc_design.capacity.fix(dfc_size[ps])

        # Solve the model
        sol = solver.solve(m, tee=True)

        # Write results to file
        _filename = "validation/" + price_signal_list[ps]["dataset"]
        _write_results(m, sol, filename=_filename) 

    return


if __name__ == "__main__":
    # single_case()
    validation_runs()