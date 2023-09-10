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
Contains functions that implement startup and shutdown constraints for the power cycle.

In this version, we implement the following and startup and shutdown procedure.

Cycle startup:
t = 0 - 1 => Ramp up fuel from 0 to 5.8518 kg/s. Power production = 0 MW
t = 1 - 2 => Ramp up fuel from 5.8515 to 10.36 kg/s. Power production = 0 MW
t = 2 - 3 => Cycle operates at P_min, (20% of the maximum load)
t = 3 - 4 => Normal operation of the plant. Cycle could ramp up to P_max

Cycle shutdown:
t = 0 - 1 => Cycle needs to ramp down to P_min
t = 1 - 2 => Cycle operates at P_min, (20% of the maximum load)
t = 2 - 3 => Ramp down fuel from 10.36 to 5.85 kg/s. Power production = 0 MW
t = 3 - 4 => Ramp down fuel from 5.8518 kg/s to 0. Power production = 0 MW
"""

from pyomo.environ import Constraint


def dfc_startup_shutdown_constraints(m):
    # The startup and shutdown constraints will be added in a separate block 
    # for eCach of the processes.

    # Create an alias for the parent block. 
    pb = m.parent_block()               # This is the multiperiod model object
    set_period = pb.set_period
    num_time_periods = len(set_period)
    des_mdl = pb.parent_block()         # This is the pyomo model object
    _params = des_mdl.dfc_design.config.model_params
    DFC_OFFLOAD = _params["op_capacity_range"][0]

    # Relationship between binary variables. This constraint implies the following
    # constraints:
    # 1. Cycle can operate at time t, only if one of the following two is true
    #     - Cycle was in operation at time t - 1
    #     - Cycle startup was initiated at time t - 2
    # 2. If the plant was in operation at time t-1 and *not* in operation at t, then plant 
    # shutdown must have been initiated at time t
    # These constraints enable us to relax binary requirement on shutdown variable
    @m.Constraint(set_period)
    def binary_var_relation(blk, t):
        if t == 1:
            # At t = 1, the generator is either running, or startup is initiated
            # Avoid shutdown at t = 1
            return pb.period[t].fs.dfc.shutdown == 0
        
        elif t == 2:
            return (
                pb.period[t].fs.dfc.op_mode - pb.period[t - 1].fs.dfc.op_mode
                == - pb.period[t].fs.dfc.shutdown
            )
        
        else:
            return (
                pb.period[t].fs.dfc.op_mode - pb.period[t - 1].fs.dfc.op_mode
                == pb.period[t - 2].fs.dfc.startup - pb.period[t].fs.dfc.shutdown
            )
    
    # If cycle startup or shutdown is initiated at time t, then no action must be
    # taken at t + 1
    @m.Constraint(set_period)
    def startup_shutdown_con_1(blk, t):
        if t + 1 > num_time_periods:
            return Constraint.Skip

        return (
            pb.period[t + 1].fs.dfc.startup + pb.period[t + 1].fs.dfc.op_mode + 
            pb.period[t + 1].fs.dfc.shutdown <= des_mdl.dfc_design.build_dfc 
            - pb.period[t].fs.dfc.startup - pb.period[t].fs.dfc.shutdown
        )

    # If cycle startup is initiated at time t i.e., startup[t] = 1, then
    # at t + 2, op_mode must be 1 and power output must be P_min
    # at t + 3, op_mode must be 1
    @m.Constraint(set_period)
    def startup_con_2(blk, t):
        if t + 2 > num_time_periods:
            return Constraint.Skip

        return pb.period[t].fs.dfc.startup <= pb.period[t + 2].fs.dfc.op_mode

    @m.Constraint(set_period)
    def startup_con_3(blk, t):
        if t + 3 > num_time_periods:
            return Constraint.Skip

        return pb.period[t].fs.dfc.startup <= pb.period[t + 3].fs.dfc.op_mode

    # Ensure that the power output is P_min at t + 2
    @m.Constraint(set_period)
    def startup_con_4(blk, t):
        if t + 2 > num_time_periods:
            return Constraint.Skip

        return (
            pb.period[t + 2].fs.dfc.power - DFC_OFFLOAD * des_mdl.dfc_design.capacity <=
            (1 - DFC_OFFLOAD) * (des_mdl.dfc_design.capacity - pb.period[t].fs.dfc.startup_capacity)
        )

    # If cycle shutdown is initiated at time t, i.e., shutdown[t] = 1, then
    # at time t - 2, op_mode must be 1 and power produced must be P_min
    # at time t - 1, op_mode must be 1 and power produced must be P_min
    # at time t, shutdown is initiated so power must be zero
    @m.Constraint(set_period)
    def shutdown_con_1(blk, t):
        if t <= 2:
            return Constraint.Skip

        return pb.period[t - 2].fs.dfc.op_mode >= pb.period[t].fs.dfc.shutdown

    @m.Constraint(set_period)
    def shutdown_con_2(blk, t):
        if t <= 1:
            return Constraint.Skip

        return pb.period[t - 1].fs.dfc.op_mode >= pb.period[t].fs.dfc.shutdown

    # Ensure that the power output is P_min at t - 2
    @m.Constraint(set_period)
    def shutdown_con_4(blk, t):
        if t <= 2:
            return Constraint.Skip

        return (
            pb.period[t - 2].fs.dfc.power - DFC_OFFLOAD * des_mdl.dfc_design.capacity <=
            (1 - DFC_OFFLOAD) * (des_mdl.dfc_design.capacity - pb.period[t].fs.dfc.shutdown_capacity)
        )

    # Ensure that the power output is P_min at t - 1
    @m.Constraint(set_period)
    def shutdown_con_5(blk, t):
        if t <= 1:
            return Constraint.Skip

        return (
            pb.period[t - 1].fs.dfc.power - DFC_OFFLOAD * des_mdl.dfc_design.capacity <=
            (1 - DFC_OFFLOAD) * (des_mdl.dfc_design.capacity - pb.period[t].fs.dfc.shutdown_capacity)
        )
        
    _min_fuel = ((_params["op_capacity_range"][0] - _params["op_curve_coeff"][1]) * _params["ng_flow_coeff"]
                 / _params["op_curve_coeff"][0])
    @m.Constraint(set_period)
    def su_sd_fuel_req(blk, t):
        # Minimum fuel requirement calculation
        k1 = 0.25 * _min_fuel * pb.period[t].fs.dfc.startup_capacity 
        k2 = 0.75 * _min_fuel * pb.period[t - 1].fs.dfc.startup_capacity if t-1 > 0 else 0
        k3 = 0.75 * _min_fuel * pb.period[t].fs.dfc.shutdown_capacity
        k4 = 0.25 * _min_fuel * pb.period[t - 1].fs.dfc.shutdown_capacity if t-1 > 0 else 0

        return pb.period[t].fs.dfc.su_sd_ng_flow == k1 + k2 + k3 + k4
