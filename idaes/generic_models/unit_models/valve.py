##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
This provides standard valve models for adaibatice control valves.  Beyond the
most common valve models, and adiabatic valve model can be added by supplying
custom callbacks for the pressure-flow relation or valve function.
"""

__Author__ = "John Eslick"

from enum import Enum

import pyomo.environ as pyo
from pyomo.common.config import ConfigValue, In

from idaes.core import declare_process_block_class
from idaes.generic_models.unit_models.pressure_changer import (
    PressureChangerData,
    ThermodynamicAssumption,
    MaterialBalanceType,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

_log = idaeslog.getLogger(__name__)


class ValveFunctionType(Enum):
    linear = 1
    quick_opening = 2
    equal_percentage = 3


def linear_cb(b):
    """
    Linear opening valve function callback.
    """
    @b.Expression(b.flowsheet().config.time)
    def valve_function(b2, t):
        return b2.valve_opening[t]


def quick_cb(b, t):
    """
    Quick opening valve function callback.
    """
    @b.Expression(b.flowsheet().config.time)
    def valve_function(b2, t):
        return sqrt(b2.valve_opening[t])


def equal_percentage_cb(b):
    """
    Equal percentage valve function callback.
    """
    b.alpha = pyo.Var(initialize=100, doc="Valve function parameter")
    b.alpha.fix()
    @b.Expression(b.flowsheet().config.time)
    def valve_function(b2, t):
        return b2.alpha ** (b2.valve_opening[t] - 1)


def pressure_flow_default_callback(b):
    """
    Add the default pressure flow relation constraint.  This will be used in the
    valve model, a custom callback is provided.
    """
    umeta = b.config.property_package.get_metadata().get_derived_units

    b.Cv = pyo.Var(
        initialize=0.1,
        doc="Valve flow coefficent",
        units=umeta("amount")/umeta("time")/umeta("pressure")**0.5
    )
    b.Cv.fix()

    b.flow_var = pyo.Reference(b.control_volume.properties_in[:].flow_mol)
    b.pressure_flow_equation_scale = lambda x : x**2

    @b.Constraint(b.flowsheet().config.time)
    def pressure_flow_equation(b2, t):
        Po = b2.control_volume.properties_out[t].pressure
        Pi = b2.control_volume.properties_in[t].pressure
        F = b2.control_volume.properties_in[t].flow_mol
        Cv = b2.Cv
        fun = b2.valve_function[t]
        return F ** 2 == Cv ** 2 * (Pi - Po) * fun ** 2


def _define_config(config):
    config.compressor = False
    config.get("compressor")._default = False
    config.get("compressor")._domain = In([False])
    config.material_balance_type = MaterialBalanceType.componentTotal
    config.get("material_balance_type")._default = \
        MaterialBalanceType.componentTotal
    config.thermodynamic_assumption = ThermodynamicAssumption.adiabatic
    config.get("thermodynamic_assumption")._default = \
        ThermodynamicAssumption.adiabatic
    config.get("thermodynamic_assumption")._domain = In(
        [ThermodynamicAssumption.adiabatic]
    )
    config.declare(
        "valve_function_callback",
        ConfigValue(
            default=ValveFunctionType.linear,
            description="Valve function type or callback for custom",
            doc="""This takes either an enumerated valve function type in: {
ValveFunctionType.linear, ValveFunctionType.quick_opening,
ValveFunctionType.equal_percentage, ValveFunctionType.custom} or a callback
function that takes a valve model object as an argument and adds a time-indexed
valve_function expression to it. Any addtional required variables, expressions,
or constraints required can also be added by the callback.""",
        ),
    )
    config.declare(
        "pressure_flow_callback",
        ConfigValue(
            default=pressure_flow_default_callback,
            description="Callback function providing the valve_function expression",
            doc="""This callback function takes a valve model object as an argument
an adds a time-indexed valve_function expression to it.  Any addtional required
variables, expressions, or constraints required can also be added by the callback.""",
        ),
    )


@declare_process_block_class("Valve", doc="Adiabatic valves")
class ValveData(PressureChangerData):
    # Same settings as the default pressure changer, but force to expander with
    # isentropic efficiency
    CONFIG = PressureChangerData.CONFIG()
    _define_config(CONFIG)

    def build(self):
        super().build()

        self.valve_opening = pyo.Var(
            self.flowsheet().config.time,
            initialize=1,
            doc="Fraction open for valve from 0 to 1",
        )
        self.valve_opening.fix()
        # If the valve function callback is set to one of the known enumerated
        # types, set the callback appropriately.  If not callable and not a known
        # type raise ConfigurationError.
        vfcb = self.config.valve_function_callback
        if not callable(vfcb):
            if vfcb == ValveFunctionType.linear:
                self.config.valve_function_callback = linear_cb
            elif vfcb == ValveFunctionType.quick_opening:
                self.config.valve_function_callback = quick_cb
            elif vfcb == ValveFunctionType.equal_percentage:
                self.config.valve_function_callback = equal_percentage_cb
            else:
                raise ConfigurationError("Invalid valve function callback.")
        self.config.valve_function_callback(self)
        self.config.pressure_flow_callback(self)

    def initialize(
        self,
        state_args={},
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg={"tol": 1e-6, "max_iter": 30},
    ):
        """
        Initialize the valve based on a deltaP guess.

        Args:
            state_args (dict): Initial state for property initialization
            outlvl : sets output level of initialization routine
            solver (str): Solver to use for initialization
            optarg (dict): Solver arguments dictionary
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")

        # storage settings what's fixed/free what's active/inactive and values
        # only for orginally fixed things.
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        self.deltaP[:].unfix()
        self.ratioP[:].unfix()

        # fix inlet and free outlet
        for t in self.flowsheet().config.time:
            for k, v in self.inlet.vars.items():
                v[t].fix()
            for k, v in self.outlet.vars.items():
                v[t].unfix()
            # to calculate outlet pressure
            Pout = self.outlet.pressure[t]
            Pin = self.inlet.pressure[t]
            if self.deltaP[t].value is not None:
                prdp = pyo.value((self.deltaP[t] - Pin) / Pin)
            else:
                prdp = -100  # crazy number to say don't use deltaP as guess
            if pyo.value(Pout / Pin) > 1 or pyo.value(Pout / Pin) < 0.0:
                if pyo.value(self.ratioP[t]) <= 1 and pyo.value(self.ratioP[t]) >= 0:
                    Pout.value = pyo.value(Pin * self.ratioP[t])
                elif prdp <= 1 and prdp >= 0:
                    Pout.value = pyo.value(prdp * Pin)
                else:
                    Pout.value = pyo.value(Pin * 0.95)
            self.deltaP[t] = pyo.value(Pout - Pin)
            self.ratioP[t] = pyo.value(Pout / Pin)

        # Make sure the initialization problem has no degrees of freedom
        # This shouldn't happen here unless there is a bug in this
        dof = degrees_of_freedom(self)
        try:
            assert dof == 0
        except:
            init_log.exception("degrees_of_freedom = {}".format(dof))
            raise

        # one bad thing about reusing this is that the log messages aren't
        # really compatible with being nested inside another initialization
        super().initialize(
            state_args=state_args, outlvl=outlvl, solver=solver, optarg=optarg
        )

        # reload original spec
        from_json(self, sd=istate, wts=sp)

    def calculate_scaling_factors(self):
        """
        Calculte pressure flow constraint scaling from flow variable scale.
        """
        # The valve of the valve opening and the output of the valve function
        # expression are between 0 and 1, so the only thing that needs to be
        # scaled here is the pressure-flow constraint, which can be scaled by
        # using the flow varable scale.  The flow variable could be defined
        # in differnt ways, so the flow variable is determined here from a
        # "flow_var[t]" reference set in the pressure-flow callback. The flow
        # term could be in verious forms, so an optional
        # "pressure_flow_equation_scale" function can be defined in the callback
        # as well.  The pressure flow function could be flow = f(Pin, Pout), but
        # it could also be flow**2 = f(Pin, Pout), ... The so
        # "pressure_flow_equation_scale" provides the form of the LHS side as
        # a function of the flow variable.

        super().calculate_scaling_factors()

        # Do some error trapping.
        if not hasattr(self, "pressure_flow_equation"):
            raise AttributeError(
                "Pressure-flow callback must define pressure_flow_equation[t]")
        # Check for flow term form if none assume flow = f(Pin, Pout)
        if hasattr(self, "pressure_flow_equation_scale"):
            ff = self.pressure_flow_equation_scale
        else:
            ff = lambda x : x
        # if the "flow_var" is not set raise an exception
        if not hasattr(self, "flow_var"):
            raise AttributeError(
                "Pressure-flow callback must define flow_var[t] reference")

        # Calculate and set the pressure-flow relation scale.
        if hasattr(self, "pressure_flow_equation"):
            for t, c in self.pressure_flow_equation.items():
                iscale.set_scaling_factor(
                    c,
                    ff(iscale.get_scaling_factor(
                        self.flow_var[t],
                        default=1,
                        warning=True)))

    def _get_performance_contents(self, time_point=0):
        pc = super()._get_performance_contents(time_point=time_point)

        pc["vars"]["Opening"] = self.valve_opening[time_point]
        try:
            pc["vars"]["Valve Coefficient"] = self.Cv
        except AttributeError:
            pass
        if self.config.valve_function == ValveFunctionType.equal_percentage:
            pc["vars"]["alpha"] = self.alpha
        return pc
