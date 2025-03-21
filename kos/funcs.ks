CLEARSCREEN.

// ############## Parameters ##############
global nd to nextnode.
global g0 to constant():G0. // gravity constant
global stage_mass to 0.607. // mass of the spin burn stage


// ############## Functions ##############
function print_status {
    parameter info.
    print "status: " +info at (5, 5).
}

function print_burn_time {
    parameter t_burn.
    print " Burn time: " + ROUND(-1*t_burn,1) + " s      " at (0,7).
}

function print_alignment {
    local align to vang(nd:deltav, ship:facing:vector).
    if align < 0.25 print "         Alignment OK: " + ROUND(align,2) + "deg    " AT(0,9).
    else print "            Alignment: " + ROUND(align,2) + "deg    " AT(0,9).
}

local active_engines is list().
function check_fuel_stability{
    local ready to 1.
    for engine in active_engines {
        set ready to ready * engine:fuelstability.
    }
    return ready.
}

function calc_spool_time{
    parameter engine.
    if engine:HasModule("ModuleEnginesRF")
    {
        local engMod is engine:GetModule("ModuleEnginesRF").
        if engMod:HasField("effective spool-up time")
            return engMod:Getfield("effective spool-up time").
    }
    return 0.01.
}

global function calc_burn_time {
    parameter deltav.  // delta-v to burn (ex: half of total_dv)
    local thrust to 0.
    local isp to 0.
    for engine in active_engines {
        set thrust to thrust + engine:maxthrust.
        set isp to isp + engine:visp * engine:maxthrust. // isp weighted by thrust
    }
    set isp to isp / thrust. // average isp by thrust

    // time = (m0 * Isp * g0 / thrust) * [1 - exp(-deltav / (Isp*g0))]
    set timeToBurn to (stage_mass * isp * g0 / thrust) * (1 - constant():e^(-deltav / (isp*g0))).
    return timeToBurn.
}

list engines in engine_list.
active_engines:add(engine_list[0]).
print active_engines[0].

print calc_burn_time(3130.6/2).
