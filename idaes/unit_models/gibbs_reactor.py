##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Standard IDAES Gibbs reactor model.
"""
from __future__ import division

# Import Pyomo libraries
from pyomo.environ import log, Reals,  Var
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume0D,
                        declare_process_block_class,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block

__author__ = "Jinliang Ma, Andrew Lee"


@declare_process_block_class("GibbsReactor")
class GibbsReactorData(UnitBlockData):
    """
    Standard Gibbs Reactor Unit Model Class

    This model assume all possible reactions reach equilibrium such that the
    system partial molar Gibbs free energy is minimized.
    Since some species mole flow rate might be very small,
    the natural log of the species molar flow rate is used.
    Instead of specifying the system Gibbs free energy as an objective
    function, the equations for zero partial derivatives of the grand function
    with Lagrangian multiple terms with repect to product species mole flow
    rates and the multiples are specified as constraints.
    """
    CONFIG = ConfigBlock()
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Gibbs reactors do not support dynamic models, thus this must be
False."""))
    CONFIG.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.enthalpyTotal,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single ethalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - ethalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))
    CONFIG.declare("momentum_balance_type", ConfigValue(
        default=MomentumBalanceType.pressureTotal,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}"""))
    CONFIG.declare("has_heat_transfer", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Heat transfer term construction flag",
        doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}"""))
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(GibbsReactorData, self).build()

        # Build Control Volume
        self.control_volume = ControlVolume0D(default={
                "dynamic": False,
                "property_package": self.config.property_package,
                "property_package_args": self.config.property_package_args})

        self.control_volume.add_state_blocks()

        self.control_volume.add_total_element_balances(
            dynamic=False,
            has_holdup=False)

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            dynamic=False,
            has_holdup=False,
            has_heat_transfer=self.config.has_heat_transfer,
            has_work_transfer=False)

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            dynamic=self.config.dynamic,
            has_holdup=False,
            has_pressure_change=self.config.has_pressure_change)

        # Add performance equations
        self.add_performance_equations()

        # Add Ports
        self.add_inlet_port()
        self.add_outlet_port()

    def add_performance_equations(self):
        """
        Define constraints which describe the behaviour of the unit model.

        Args:
            None

        Returns:
            None
        """
        object.__setattr__(self,
                           "component_list",
                           self.control_volume.component_list)
        object.__setattr__(self,
                           "phase_list",
                           self.control_volume.phase_list)
        object.__setattr__(self,
                           "element_list",
                           self.control_volume.element_list)

        # Add Lagrangian multiplier variables
        self.lagrange_mult = Var(self.time,
                                 self.element_list,
                                 domain=Reals,
                                 initialize=1000,
                                 doc="Lagrangian multipliers")

        # Use Lagrangian multiple method to derive equations for Out_Fi
        # Use RT*lagrange as the Lagrangian multiple such that lagrange is in
        # a similar order of magnitude as log(Yi)

        @self.Constraint(self.time,
                         self.phase_list,
                         self.component_list,
                         doc="Gibbs energy minimisation constraint")
        def gibbs_minimization(b, t, p, j):
            # Use natural log of species mole flow to avoid Pyomo solver
            # warnings of reaching infeasible point
            return 0 == (
                b.control_volume.properties_out[t].gibbs_mol_phase_comp[p, j] +
                b.control_volume.properties_out[t].gas_const *
                b.control_volume.properties_out[t].temperature *
                (log(b.control_volume.properties_out[t].mole_frac[j]) +
                 sum(b.lagrange_mult[t, e] *
                     b.control_volume.properties_out[t].element_comp[j][e]
                     for e in b.element_list)))

        # Set references to balance terms at unit level
        if (self.config.has_heat_transfer is True and
                self.config.energy_balance_type != EnergyBalanceType.none):
            object.__setattr__(self,
                               "heat",
                               self.control_volume.heat)
