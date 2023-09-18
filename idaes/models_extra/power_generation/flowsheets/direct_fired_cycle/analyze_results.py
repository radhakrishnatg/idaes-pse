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

import json
import pandas as pd
import matplotlib.pyplot as plt
from pyomo.environ import value


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

    11: {"dataset": "NETL_ERCOT_GEN_0", "carbon_tax": 0},
    12: {"dataset": "NETL_ERCOT_GEN_25", "carbon_tax": 25},
    13: {"dataset": "NETL_ERCOT_GEN_50", "carbon_tax": 50},
    14: {"dataset": "NETL_ERCOT_GEN_100", "carbon_tax": 100},
    15: {"dataset": "NETL_ERCOT_GEN_250", "carbon_tax": 250},

    16: {"dataset": "ISO_NE_2022", "carbon_tax": 0},
    17: {"dataset": "MISO_INDIANA_2022", "carbon_tax": 0},
    18: {"dataset": "SPP_NORTH_2022", "carbon_tax": 0},
    19: {"dataset": "SPP_SOUTH_2022", "carbon_tax": 0},

    "V1": {"dataset": "EXEMPLAR_LARGE", "carbon_tax": 0},
    "V2": {"dataset": "EXEMPLAR_MEDIUM", "carbon_tax": 0},
    "V3": {"dataset": "EXEMPLAR_SMALL", "carbon_tax": 0},
}


def _write_results(
    m,
    sol,
    lox_withdrawal=False,
    includes_nlu=False,
    filename="results",
):
    _filename = filename + ".xlsx"

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
        results["Power_dfc_to_tank"] = [value(fs.power_dfc_to_tank) for fs in set_flowsheets]
        results["Power_grid_to_nlu"] = [value(fs.power_grid_to_nlu) for fs in set_flowsheets]
        results["Power_grid_to_tank"] = [value(fs.power_grid_to_tank) for fs in set_flowsheets]

    if lox_withdrawal:
        results["GOx_fraction"] = [value(fs.gox_fraction) for fs in set_flowsheets]

        results["Tank_init_holdup"] = [value(fs.tank.initial_holdup) for fs in set_flowsheets]
        results["Tank_final_holdup"] = [value(fs.tank.final_holdup) for fs in set_flowsheets]
        results["Tank_lox_in"] = [value(fs.tank.lox_in) for fs in set_flowsheets]
        results["Tank_lox_out"] = [value(fs.tank.lox_out) for fs in set_flowsheets]

        results["Power_dfc_to_tank"] = [value(fs.power_dfc_to_tank) for fs in set_flowsheets]
        results["Power_grid_to_tank"] = [value(fs.power_grid_to_tank) for fs in set_flowsheets]

    results_df = pd.DataFrame(results)
    results_df.to_excel(_filename)

    lower_bnd = sol["Problem"][0]["Lower bound"]
    upper_bnd = sol["Problem"][0]["Upper bound"]

    elec_rev = [value(fs.electricity_revenue) for fs in set_flowsheets]
    elec_cost = [value(fs.electricity_cost) for fs in set_flowsheets]
    fuel_cost = [value(fs.fuel_cost) for fs in set_flowsheets]
    co2_price = [value(fs.co2_price) for fs in set_flowsheets]
    dfc_vom = [value(fs.dfc.non_fuel_vom) for fs in set_flowsheets]
    asu_vom = [value(fs.asu.non_fuel_vom) for fs in set_flowsheets]
    if includes_nlu:
        nlu_vom = [value(fs.nlu.non_fuel_vom) for fs in set_flowsheets]
    else:
        nlu_vom = [0]

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
        "DFC CAPEX": value(m.dfc_design.capex) / 1e3,
        "ASU CAPEX": value(m.asu_design.capex) / 1e3,
        "NLU CAPEX": value(m.nlu_design.capex) / 1e3 if hasattr(m, "nlu_Design") else 0,
        "TANK CAPEX": value(m.tank_design.capex) / 1e3 if hasattr(m, "tank_design") else 0,
        "Total CAPEX": m.CAPEX.value / 1e3,
        "Total FOM": m.FOM.value / 1e3,
        "Total Revenue": sum(elec_rev) / 1e3,
        "Total Elec Cost": sum(elec_cost) / 1e3,
        "Total Fuel Cost": sum(fuel_cost) / 1e3,
        "Total CO2 Price": sum(co2_price) / 1e3,
        "Total VOM": (sum(dfc_vom) + sum(asu_vom) + sum(nlu_vom)) / 1e3,
        "Total Corp. tax": m.CORP_TAX.value / 1e3,
        "Total Net Profit": m.NET_PROFIT.value / 1e3,

        "Num Vars": sol["Problem"][0]["Number of variables"],
        "Num Bin Vars": sol["Problem"][0]["Number of binary variables"],
        "Num constraints": sol["Problem"][0]["Number of constraints"],
    }

    with open(filename + ".json", "w") as fp:
        json.dump(solution, fp, indent=4)


def summarize_sa_results(folder, res_filename="results.xlsx"):
    _price_signal_list = [price_signal_list[i]["dataset"] for i in range(1, 20)]
    
    key_list = [
        "NPV", "Upper bound", "Gap", "DFC_Capacity", "ASU_Capacity", "NLU_Capacity", "Tank_Capacity",
        "DFC Startups", "DFC Shutdowns", "DFC Capacity Factor", "ASU Startups", "ASU Shutdowns", 
        "DFC CAPEX", "ASU CAPEX", "NLU CAPEX", "TANK CAPEX", "Total CAPEX", "Total FOM", "Total Revenue", 
        "Total Elec Cost", "Total Fuel Cost", "Total CO2 Price", "Total VOM", "Total Corp. tax", 
        "Total Net Profit", "Num Vars", "Num Bin Vars", "Num constraints",
    ]

    results_summary = {}
    for k in key_list:
        results_summary[k] = []  

    for ps in _price_signal_list:
        filename = folder + ps + ".json"

        with open(filename) as fp:
            data = json.load(fp)

        for k in key_list:
            if k == "Upper bound":
                results_summary[k].append(data[k] / 1000)

            elif k == "NLU_Capacity" or k == "Tank_Capacity":
                results_summary[k].append(0 if data[k] == "N/A" else data[k])

            elif k == "DFC Capacity Factor":
                results_summary[k].append(data["Total Power Produced"] / (8760 * data["DFC_Capacity"]))

            else:
                results_summary[k].append(data[k])

    results_df = pd.DataFrame(results_summary, index=_price_signal_list)
    results_df.to_excel(res_filename)
