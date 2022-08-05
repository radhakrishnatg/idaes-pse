from pyomo.environ import (
    Var,
    NonNegativeReals,
    Constraint,
)
from idaes.core import FlowsheetBlock
from .unit_models import (
    DFCOperation,
    MonoASUOperation,
    NLUOperation,
    OxygenTankOperation,
)


def build_dfc_flowsheet(
    m,
    dfc_design,
    asu_design,
    nlu_design,
    tank_design
):
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Append unit models
    m.fs.dfc = DFCOperation(default={"design_blk": dfc_design})
    m.fs.asu = MonoASUOperation(default={"design_blk": asu_design})
    m.fs.nlu = NLUOperation(default={"design_blk": nlu_design})
    m.fs.tank = OxygenTankOperation(default={"design_blk": tank_design})

    # Declare new variables
    m.fs.power_dfc_to_grid = Var(within=NonNegativeReals)
    m.fs.power_dfc_to_asu = Var(within=NonNegativeReals)
    m.fs.power_dfc_to_nlu = Var(within=NonNegativeReals)
    m.fs.power_dfc_to_tank = Var(within=NonNegativeReals)
    m.fs.power_grid_to_asu = Var(within=NonNegativeReals)
    m.fs.power_gird_to_nlu = Var(within=NonNegativeReals)
    m.fs.power_grid_to_tank = Var(within=NonNegativeReals)

    m.fs.oxygen_asu_to_dfc = Var(within=NonNegativeReals)

    # Power balance across the power cycle
    m.fs.dfc_power_balance = Constraint(
        expr=m.fs.dfc.power == m.fs.power_dfc_to_grid + m.fs.power_dfc_to_asu + 
        m.fs.power_dfc_to_nlu + m.fs.power_dfc_to_tank,
    )

    # Power balance across the ASU
    m.fs.asu_power_balance = Constraint(
        expr=m.fs.asu.power == m.fs.power_dfc_to_asu + m.fs.power_grid_to_asu
    )

    # Power balance across the NLU
    m.fs.nlu_power_balance = Constraint(
        expr=m.fs.nlu.power == m.fs.power_dfc_to_nlu + m.fs.power_gird_to_nlu
    )

    # Power balance across the tank
    m.fs.tank_power_balance = Constraint(
        expr=m.fs.tank.power == m.fs.power_dfc_to_tank + m.fs.power_grid_to_tank
    )

    # Oxygen balance across the power cycle
    m.fs.dfc_oxygen_balance = Constraint(
        expr=m.fs.dfc.o2_flow == m.fs.oxygen_asu_to_dfc + m.fs.tank.Lox_out
    )

    # Oxygen balance across the ASU
    m.fs.asu_oxygen_balance = Constraint(
        expr=m.fs.asu.o2_flow == m.fs.oxygen_asu_to_dfc + m.fs.nlu.o2_flow
    )
