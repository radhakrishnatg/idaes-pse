"""
Contains functions that implement startup and shutdown constraints for the ASU

In this version, we implement the following and startup and shutdown procedure.

Cycle startup:
t = 0 - 8 => Consumes 80% of the max_power, but does not produce any oxygen
t > 8 => Can operate anywhere between P_min and P_max

Cycle shutdown:
t = 0 - 1 => Consumes 50% of the max_power and does not produce any oxygen
"""

from pyomo.environ import Constraint

# Startup time for the ASU
SU_TIME = 8


def asu_startup_shutdown_constraints(m):
    # The startup and shutdown constraints will be added in a separate block 
    # for eCach of the processes.

    # Create an alias for the parent block. 
    pb = m.parent_block()               # This is the multiperiod model object
    set_period = pb.set_period
    num_time_periods = len(set_period)
    des_mdl = pb.parent_block()         # This is the pyomo model object
    _power_req_coeff = des_mdl.asu_design.config.model_params["power_req_coeff"]

    # Relationship between binary variables. This constraint implies the following
    # constraints:
    # 1. Cycle can operate at time t, only if one of the following two is true
    #     - Cycle was in operation at time t - 1
    #     - Cycle startup was initiated at time t - 8
    # 2. If the plant was in operation at time t-1 and *not* in operation at t, then plant 
    # shutdown must have been initiated at time t
    # These constraints enable us to relax binary requirement on shutdown variable
    @m.Constraint(set_period)
    def binary_var_relation(blk, t):
        if t == 1:
            # At t = 1, the generator is either running, or startup is initiated
            # Avoid shutdown at t = 1
            return pb.period[t].fs.asu.shutdown == 0
        
        else:
            return (
                pb.period[t].fs.asu.op_mode - pb.period[t - 1].fs.asu.op_mode
                == (pb.period[t - SU_TIME].fs.asu.startup if t-SU_TIME > 0 else 0)
                - pb.period[t].fs.asu.shutdown
            )

    # After plant startup, the unit cannot produce oxygen in the next 8 hours, but
    # consumes 80% of P_max power during that period
    @m.Constraint(set_period, [j for j in range(1, SU_TIME)])
    def startup_con_1(blk, t, i):
        if t + i > num_time_periods:
            return Constraint.Skip

        return (
            pb.period[t + i].fs.asu.startup + pb.period[t + i].fs.asu.op_mode +
            pb.period[t + i].fs.asu.shutdown <= 1 - pb.period[t].fs.asu.startup
        )

    # I think this constraint will automatically be satisfied by the optimal solution
    # so, we do not have to add it explicitly. 
    # ASU must operate at the ninth hour after startup is initiated
    @m.Constraint(set_period)
    def startup_con_3(blk, t):
        if t + SU_TIME > num_time_periods:
            return Constraint.Skip

        return pb.period[t + SU_TIME].fs.asu.op_mode >= pb.period[t].fs.asu.startup

    # Shutdown can be initiated at time t only if the ASU was operating at time t - 1
    @m.Constraint(set_period)
    def shutdown_con_1(blk, t):
        if t == 1:
            return Constraint.Skip

        # I suspect this will always be satisfied by the optimal solution. Need to verify it though
        return pb.period[t].fs.asu.shutdown <= pb.period[t - 1].fs.asu.op_mode

    # During the shutdown, power consumption is 50% of the P_max and O2 is not produced
    @m.Constraint(set_period)
    def startup_shutdown_power(blk, t):
        return (
            pb.period[t].fs.asu.su_sd_power == 
            0.5 * _power_req_coeff * pb.period[t].fs.asu.shutdown_o2_flow +
            0.8 * _power_req_coeff * sum(
                pb.period[t - i].fs.asu.startup_o2_flow 
                for i in range(SU_TIME) if t - i > 0
            )
        )
