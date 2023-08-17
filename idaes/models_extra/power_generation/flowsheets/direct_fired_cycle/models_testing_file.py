from pyomo.environ import ConcreteModel
from idaes.apps.grid_integration import MultiPeriodModel
from idaes.power_generation.flowsheets.direct_fired_cycle.unit_models import (
    get_lmp_data,
    DFCDesign,
    MonoASUDesign,
    NLUDesign,
    OxygenTankDesign,
)
from idaes.power_generation.flowsheets.direct_fired_cycle.dfc_flowsheet import build_dfc_flowsheet


def get_linking_var_pairs(m1, m2):
    return [(m1.fs.tank.final_holdup, m2.fs.tank.initial_holdup)]


def design_model_testing():
    m = ConcreteModel()
    get_lmp_data(m, dataset="NREL", location="PJM-W", carbon_tax=150)

    m.dfc_design = DFCDesign()
    m.asu_design = MonoASUDesign()
    m.nlu_design = NLUDesign()
    m.tank_design = OxygenTankDesign()

    
    # Currently, we are not varying the capacity of the power cycle, so fixing the value
    # of capacity. Also, we want the power cycle to be built, so fixing build_dfc. 
    m.dfc_design.build_dfc.fix(1)
    m.dfc_design.capacity.fix()

    m.mp_model = MultiPeriodModel(
        n_time_points=10,
        process_model_func=build_dfc_flowsheet,
        linking_variable_func=get_linking_var_pairs,
        use_stochastic_build=True,
        flowsheet_options={
            "dfc_design": m.dfc_design,
            "asu_design": m.asu_design,
            "nlu_design": m.nlu_design,
            "tank_design": m.tank_design,
        },
    )

    # build_dfc_flowsheet(
    #     m,
    #     dfc_design=m.dfc_design,
    #     asu_design=m.asu_design,
    #     nlu_design=m.nlu_design,
    #     tank_design=m.tank_design,
    # )

    return m
