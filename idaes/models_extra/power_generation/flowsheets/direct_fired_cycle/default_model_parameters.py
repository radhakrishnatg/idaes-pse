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
This script contains parameter data

The variable O&M cost (non-fuel) is calculated as follows:
Non-fuel VOM includes: TBD.
Total VOM = 3339.92 + 18048.07 = 21387.99. Split the VOM for inidividual units
in the ratio of their CAPEX. 1128855 / (1128855 + 545522) = 0.6742.

The VOM of the DFC cycle = 0.6742 * 21387.99 = 14419.78286 
 ==> 14419.78286 / (0.85 * 8760) = 1.93658 [in $1000/hr]
     This number will be scaled linearly with normalized DFC's capacity

The VOM of the ASU unit = 0.3258 * 21387.99 = 6968.207
 ==> 6968.207 / (0.85 * 8760) = 0.93583 [in $1000/hr]
     This number will be scaled linearly with normalized ASU's capacity

The VOM of the NLU unit = 5404.96
 ==> 5404.96 / (0.85 * 8760) = 0.7259 [in $1000/hr]
     This number will be scaled linearly with normalized NLU's capacity


Notes:
1. DFC's net power = Gross power - compressor power - all auxiliaries except ASU.

2. We assume that the natural gas (NG) requirement at full load is proportional to the design
   capacity. Let 'ng_flow_coeff' denote the proportinality constant. The value can be 
   fitted if multiple data points are available. But if only one data point is available,
   then it is the ratio of NG flowrate at full load and design capacity.

3. We assume that the CAPEX of DSU = k1 * design net power + k2. The net power is defined in
   Note 1. Provide k1 and k2 as a list i.e., [k1, k2]

4. The fixed O&M cost is calculated as a fraction of CAPEX, and it includes: Operating labor,
   maintenance labor, labor administration and overhead charges, and property taxes and insurance.

5. Power cycle must operate between a minimum stable load and full load. Let the minimum stable 
   load = 0.3 * full load. Then, the operating capacity range, 'op_capacity_range' = (0.3, 1)

6. For off-design performance curve relates the off-design operation. The performance curve is
   r_P = k1 * r_NG + k2, where r_P = off-design power / design power and r_NG = NG requirement
   at part load / NG requirement at full load. Provide k1 and k2 as a list i.e., [k1, k2]

7. Natural gas molar composition is assumed to be 93.1% methane, 3.2% ethane, 0.7% propane,
   0.4% butane, 1% CO2, and 1.6% N2. Therefore, one kmol of natural gas yields 1.042 kmol
   of CO2, or 1.042 * 44.01 = 45.8584 kg of CO2. Molar mass of natural gas is 17.3268 kg/kmol.
   Therefore, 1 kg of natural gas emits 45.8584/17.3268 = 2.6467 kg of CO2.

8. Non-fuel variable O&M costs include consumables (maintenance material, water, water treatment, 
   chemicals, etc.) and waste disposal. Non-fuel variable O&M is not available for individual 
   units. We split this cost proportional to the capital costs of DFC and ASU. We calculate the
   vom as the sum of const_vom and var_vom. const_vom adds a fixed cost when the power cycle is 
   operating regardless of the power output. var_vom adds a variable cost that is a function of
   total power produced. const_vom = k1 * design net power + k2 and var_vom = k1 * net power
"""

# TODO: Contact Sandeep for reference. 
# HHV of NG = 22499.17034 btu/lb = 0.0496 MMBtu/kg
NG_HHV = 0.0496
HR_TO_SEC = 3600
LBS_TO_KG = 0.453592

# Parameter data for DFC
# Calculating the FOM as 3.157% of the CAPEX. The percentage value is obtained from
# 35,641.27 (FOM) / 1,128,855 (CAPEX). The FOM is also assumed 
# to vary linearly with the capacity.
DFC_PARAMS = {
    "des_capacity_range" : (176.032, 777.685),        # [MW] Bounds on DFC's net power (See Note 1)
    "ng_flow_coeff"      : 0.034774,                  # [kg/s/MW] (See Note 2) 
    "o2_ng_ratio"        : 3.78495,                   # [-] Oxygen/NG flowratio
    "capex"              : [1067.8405, 239897.1320],  # [$1000] CAPEX of DFC (See Note 3)
    "fom_factor"         : 0.03128,                   # [-] Multiplier for FOM (See Note 4)
    "op_capacity_range"  : (0.2, 1),                  # [-] Operating range of the power cycle (See Note 5)
    "op_curve_coeff"     : [1.4811, -0.4811],         # Coefficients of performance curve (See Note 6)
    "co2_emission_rate"  : 2.6467,                    # [kg CO2/kg NG] (See Note 7).
    "co2_capture_rate"   : 0.985,                     # [-] Fraction of CO2 captured
    "var_vom_coeff"      : 0,                         # [$1000/MWh] Non-fuel VOM (See Note 8)
    "const_vom_coeff"    : [0.0018033, 0.3736672],    # [$1000/h] Non-fuel VOM (See Note 8)
}

# Parameter data for Monolithic ASU
MONO_ASU_PARAMS = {
    "des_capacity_range" : (20, 110),                 # [kg/s] Bounds on ASU capacity
    "power_req_coeff"    : 1.3086,                    # [MW] Power requirement at max capacity 
    "capex"              : [5394.1471, 17531.4670],   # [$1000] CAPEX of ASU 
    "fom_factor"         : 0.03128,                   # [-] Multiplier for FOM
    "op_capacity_range"  : (0.3, 1),                  # [-] Operating range of ASU
    "op_curve_coeff"     : [0.9625, 0.0375],          # Coefficients of performance curve
    "argon_price"        : 0,                         # [$/kg] Selling price of Argon
    "nitrogen_price"     : 0,                         # [$/kg] Selling price of nitrogen
    "var_vom_coeff"      : 0,                         # [$1000/MWh] Non-electricity VOM 
    "const_vom_coeff"    : [0.00903, 0.0205359],      # [$1000/hr] Non-electricity VOM
}

# Parameter data for ASU
ASU_PARAMS = {
    "des_capacity_range" : (5, 27.5),                 # [kg/s] Bounds on ASU capacity
    "power_req_coeff"    : 1.3086,                    # [MW] Power requirement at max capacity 
    "capex"              : [1348.536775, 4382.86675], # [$1000] CAPEX of ASU 
    "fom_factor"         : 0.03128,                   # [-] Multiplier for FOM
    "op_capacity_range"  : (0.7, 1),                  # [-] Operating range of ASU
    "op_curve_coeff"     : [0.9625, 0.0375],          # Coefficients of performance curve
    "argon_price"        : 0,                         # [$/kg] Selling price of Argon
    "nitrogen_price"     : 0,                         # [$/kg] Selling price of nitrogen
    "var_vom_coeff"      : 0,                         # [$1000/MWh] Non-electricity VOM 
    "const_vom_coeff"    : [0.0022575, 0.005133975],  # [$1000/hr] Non-electricity VOM
}

# Parameter data for NLU
NLU_PARAMS = {
    "des_capacity_range" : (20, 110),                 # [kg/s] Bounds on ASU capacity
    "power_req_coeff"    : 1.7873,                    # [MW] Power requirement at max capacity 
    "capex"              : [1856.93338, 7023.3325],   # [$1000] CAPEX of NLU 
    "fom_factor"         : 0.03128,                   # [-] Multiplier for FOM
    "op_capacity_range"  : (0.3, 1),                  # [-] Operating range of ASU
    "var_vom_coeff"      : 0,                         # [$1000/MWh] Non-electricity VOM 
    "const_vom_coeff"    : [0.005498, 0],             # [$1000/hr] Non-electricity VOM
}

# Parameter data for LOx tank
O2_TANK_PARAMS = {
    "des_capacity_range" : (2000, 400000),            # [kg/s] Bounds on ASU capacity
    "capex"              : [0.98167214, 2779.905438], # [$1000] CAPEX of tank
    "fom_factor"         : 0.03128,                   # [-] Multiplier for FOM
    "min_holdup"         : 0.1,                       # [-] Minimum holdup fraction
    "power_req_coeff"    : 1.47054 / 112.0865459,     # [MW/kg/s] Power requirement for storage
    "tank_constraint"    : "periodic",                # "Periodic" or "initial_holdup"
}


# Cashflow parameters
CASHFLOW_PARAMS = {
    "plant_life"        : 30,                         # [-] Plant lifetime in years
    "discount_rate"     : 0.075,                      # [-] Discount rate
    "tax_rate"          : 0.2,                        # [-] Corporate tax rate
    "electricity_cost"  : 5,                          # [$/MWh] Excess penalty for purchasing electricity
}
discount_rate = CASHFLOW_PARAMS["discount_rate"]
plant_life = CASHFLOW_PARAMS["plant_life"]
CASHFLOW_PARAMS["FCR"] = 1 / ((1 - (1 + discount_rate) ** (- plant_life)) / discount_rate)
