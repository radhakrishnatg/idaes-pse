##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Reboiler model for distillation.

While the reboiler model, is fairly simple, a major
portion of this code has gone into making this generic and be able to handle
different state variables and the associated splits.
"""

__author__ = "Jaffer Ghouse"

import logging

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.network import Port
from pyomo.environ import Reference, Expression, Var, Constraint, \
    TerminationCondition

# Import IDAES cores
from idaes.logger import getIdaesLogger, getInitLogger, condition
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        MaterialBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import PropertyPackageError

_log = getIdaesLogger(__name__)


@declare_process_block_class("Reboiler")
class ReboilerData(UnitModelBlockData):
    """
    Reboiler unit for distillation model.
    Unit model to reboil the liquid from the bottom tray of
    the distillation column.
    """
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare("has_boilup_ratio", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Boilup ratio term construction flag",
        doc="""Indicates whether terms for boilup ratio should be
constructed,
**default** - False.
**Valid values:** {
**True** - include construction of boilup ratio constraint,
**False** - exclude construction of boilup ratio constraint}"""))
    CONFIG.declare("material_balance_type", ConfigValue(
        default=MaterialBalanceType.useDefault,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentPhase.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}"""))
    CONFIG.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.useDefault,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
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
        """Build the model.

        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(ReboilerData, self).build()

        # Add Control Volume for the Reboiler
        self.control_volume = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.property_package,
            "property_package_args": self.config.property_package_args})

        self.control_volume.add_state_blocks(
            has_phase_equilibrium=True)

        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=True)

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=True)

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        self.boilup_ratio = Var(initialize=0.5, doc="boilup ratio for reboiler")

        if self.config.has_boilup_ratio is True:
            def rule_boilup_ratio(self, t):
                if hasattr(self.control_volume.properties_out[t],
                           "flow_mol_phase"):
                    return self.boilup_ratio * \
                        self.control_volume.properties_out[t].\
                        flow_mol_phase["Liq"] == self.control_volume.\
                        properties_out[t].flow_mol_phase["Vap"]
                elif hasattr(self.control_volume.properties_out[t],
                             "flow_mol_phase_comp"):
                    return self.boilup_ratio * \
                        sum(self.control_volume.properties_out[t].
                            flow_mol_phase_comp["Liq", i]
                            for i in self.control_volume.properties_out[t].
                            _params.component_list) == \
                        sum(self.control_volume.properties_out[t].
                            flow_mol_phase_comp["Vap", i]
                            for i in self.control_volume.properties_out[t].
                            _params.component_list)
                else:
                    raise Exception("Unsupported flow variables")
            self.eq_boilup_ratio = Constraint(self.flowsheet().time,
                                              rule=rule_boilup_ratio)

        self._make_ports()

        self._make_splits_reboiler()

        # Add object reference to variables of the control volume
        # Reference to the heat duty
        add_object_reference(self, "heat_duty", self.control_volume.heat)
        # Reference to the pressure drop (if set to True)
        if self.config.has_pressure_change:
            add_object_reference(self, "deltaP", self.control_volume.deltaP)

    def _make_ports(self):

        # Add Ports for the reboiler
        # Inlet port (the vapor from the top tray)
        self.add_inlet_port()

        # Outlet ports that always exist irrespective of reboiler type
        self.bottoms = Port(noruleinit=True, doc="Bottoms stream.")

        self.vapor_reboil = Port(noruleinit=True,
                                 doc="Vapor outlet stream that is returned to "
                                 "to the bottom tray.")

    def _make_splits_reboiler(self):
        # Get dict of Port members and names
        member_list = self.control_volume.\
            properties_out[0].define_port_members()

        # Create references and populate the reflux, distillate ports
        for k in member_list:
            # Create references and populate the intensive variables
            if "flow" not in k and "frac" not in k and "enth" not in k:
                if not member_list[k].is_indexed():
                    var = self.control_volume.properties_out[:].\
                        component(member_list[k].local_name)
                else:
                    var = self.control_volume.properties_out[:].\
                        component(member_list[k].local_name)[...]

                # add the reference and variable name to the reflux port
                self.bottoms.add(Reference(var), k)

                # add the reference and variable name to the
                # vapor outlet port
                self.vapor_reboil.add(Reference(var), k)

            elif "frac" in k and ("mole" in k or "mass" in k):

                # Mole/mass frac is typically indexed
                index_set = member_list[k].index_set()

                # if state var is not mole/mass frac by phase
                if "phase" not in k:
                    # Assuming the state block has the var
                    # "mole_frac_phase_comp". Valid if VLE is supported
                    # Create a string "mole_frac_phase_comp" or
                    # "mass_frac_phase_comp". Cannot directly append phase to
                    # k as the naming convention is phase followed by comp

                    str_split = k.split('_')
                    local_name = '_'.join(str_split[0:2]) + \
                        "_phase" + "_" + str_split[2]

                    # Rule for liquid fraction
                    def rule_liq_frac(self, t, i):
                        return self.control_volume.properties_out[t].\
                            component(local_name)["Liq", i]
                    self.e_liq_frac = Expression(
                        self.flowsheet().time, index_set,
                        rule=rule_liq_frac)

                    # Rule for vapor fraction
                    def rule_vap_frac(self, t, i):
                        return self.control_volume.properties_out[t].\
                            component(local_name)["Vap", i]
                    self.e_vap_frac = Expression(
                        self.flowsheet().time, index_set,
                        rule=rule_vap_frac)

                    # add the reference and variable name to the
                    # distillate port
                    self.bottoms.add(self.e_liq_frac, k)

                    # add the reference and variable name to the
                    # vapor port
                    self.vapor_reboil.add(self.e_vap_frac, k)

                else:

                    # Assumes mole_frac_phase or mass_frac_phase exist as
                    # state vars in the port and therefore access directly
                    # from the state block.
                    var = self.control_volume.properties_out[:].\
                        component(member_list[k].local_name)[...]

                    # add the reference and variable name to the distillate port
                    self.bottoms.add(Reference(var), k)

                    # add the reference and variable name to the boil up port
                    self.vapor_reboil.add(Reference(var), k)
            elif "flow" in k:
                if "phase" not in k:

                    # Assumes that here the var is total flow or component
                    # flow. However, need to extract the flow by phase from
                    # the state block. Expects to find the var
                    # flow_mol_phase or flow_mass_phase in the state block.

                    # Check if it is not indexed by component list and this
                    # is total flow
                    if not member_list[k].is_indexed():
                        # if state var is not flow_mol/flow_mass
                        # by phase
                        local_name = str(member_list[k].local_name) + \
                            "_phase"

                        # Rule for vap flow
                        def rule_vap_flow(self, t):
                            return self.control_volume.properties_out[t].\
                                component(local_name)["Vap"]
                        self.e_vap_flow = Expression(
                            self.flowsheet().time,
                            rule=rule_vap_flow)

                        # Rule to link the liq flow to the distillate
                        def rule_bottoms_flow(self, t):
                            return self.control_volume.properties_out[t].\
                                component(local_name)["Liq"]
                        self.e_bottoms_flow = Expression(
                            self.flowsheet().time,
                            rule=rule_bottoms_flow)

                    else:
                        # when it is flow comp indexed by component list
                        str_split = \
                            str(member_list[k].local_name).split("_")
                        if len(str_split) == 3 and str_split[-1] == "comp":
                            local_name = str_split[0] + "_" + \
                                str_split[1] + "_phase_" + "comp"

                        # Get the indexing set i.e. component list
                        index_set = member_list[k].index_set()

                        # Rule for vap flow
                        def rule_vap_flow(self, t, i):
                            return self.control_volume.properties_out[t].\
                                component(local_name)["Vap", i]
                        self.e_vap_flow = Expression(
                            self.flowsheet().time, index_set,
                            rule=rule_vap_flow)

                        # Rule to link the liq flow to the distillate
                        def rule_bottoms_flow(self, t, i):
                            return self.control_volume.properties_out[t].\
                                component(local_name)["Liq", i]
                        self.e_bottoms_flow = Expression(
                            self.flowsheet().time, index_set,
                            rule=rule_bottoms_flow)

                    # add the reference and variable name to the
                    # distillate port
                    self.bottoms.add(self.e_bottoms_flow, k)

                    # add the reference and variable name to the
                    # distillate port
                    self.vapor_reboil.add(self.e_vap_flow, k)
            elif "enth" in k:
                if "phase" not in k:
                    # assumes total mixture enthalpy (enth_mol or enth_mass)
                    if not member_list[k].is_indexed():
                        # if state var is not enth_mol/enth_mass
                        # by phase, add _phase string to extract the right
                        # value from the state block
                        local_name = str(member_list[k].local_name) + \
                            "_phase"
                    else:
                        raise PropertyPackageError(
                            "Expected an unindexed variable.")

                    # Rule for vap enthalpy. Setting the enthalpy to the
                    # enth_mol_phase['Vap'] value from the state block
                    def rule_vap_enth(self, t):
                        return self.control_volume.properties_out[t].\
                            component(local_name)["Vap"]
                    self.e_vap_enth = Expression(
                        self.flowsheet().time,
                        rule=rule_vap_enth)

                    # Rule to link the liq flow to the distillate.
                    # Setting the enthalpy to the
                    # enth_mol_phase['Liq'] value from the state block
                    def rule_bottoms_enth(self, t):
                        return self.control_volume.properties_out[t].\
                            component(local_name)["Liq"]
                    self.e_bottoms_enth = Expression(
                        self.flowsheet().time,
                        rule=rule_bottoms_enth)

                    # add the reference and variable name to the
                    # distillate port
                    self.bottoms.add(self.e_bottoms_enth, k)

                    # add the reference and variable name to the
                    # distillate port
                    self.vapor_reboil.add(self.e_vap_enth, k)
                elif "phase" in k:
                    # assumes enth_mol_phase or enth_mass_phase.
                    # This is an intensive property, you create a direct
                    # reference irrespective of the reflux, distillate and
                    # vap_outlet

                    # Rule for vap flow
                    if not member_list[k].is_indexed():
                        var = self.control_volume.properties_out[:].\
                            component(member_list[k].local_name)
                    else:
                        var = self.control_volume.properties_out[:].\
                            component(member_list[k].local_name)[...]

                    # add the reference and variable name to the distillate port
                    self.bottoms.add(Reference(var), k)

                    # add the reference and variable name to the
                    # vapor outlet port
                    self.vapor_reboil.add(Reference(var), k)
                else:
                    raise Exception(
                        "Unrecognized enthalpy state variable. "
                        "Only total mixture enthalpy or enthalpy by "
                        "phase are supported.")

    def initialize(self, solver=None, outlvl=None):

        # TODO: Fix the inlets to the reboiler to the vapor flow from
        # the top tray or take it as an argument to this method.

        init_log = getInitLogger(self.name, outlvl)

        # Initialize the inlet and outlet state blocks
        self.control_volume.initialize(outlvl=outlvl)

        if solver is not None:
            if outlvl > 2:
                tee = True
            else:
                tee = False

            solver_output = solver.solve(self, tee=tee)

            if solver_output.solver.termination_condition == \
                    TerminationCondition.optimal:
                init_log.log(4, 'Reboiler Initialisation Complete, {}.'
                             .format(condition(solver_output)))
