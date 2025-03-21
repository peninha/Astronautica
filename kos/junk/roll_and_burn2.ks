parameter total_dv, stage_mass.
parameter setup_time to 120, rcs_spin_time to 30, rcs_ullage_time to 2, stop_burn_tol to 0.4, burn_delay to 2.2.

CLEARSCREEN.
// Using Delta-V to calculate the burn

// ############## Parameters ##############
local nd to nextnode.
local g0 to constant():G0. // gravity constant

print "           burn delay: " + burn_delay at (0,18).

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

function calc_burn_time {
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

print "######### Next Node Burning Program #########" at(0, 3).
print_status("Warping to node...                        ").

// ############## Active engines calculations ##############

local ullage_time to 0.
local spool_time to 0.
list engines in engine_list.
local i to 0.
for engine in engine_list {
    if engine:ignition // Add future active engines to the list until the first currently active engine is found
        break.
    engine:activate(). // activate the engine
    active_engines:add(engine). // add the engine to the list of future active engines
	set spool_time to spool_time + calc_spool_time(engine).
    if engine:ullage {
        set ullage_time to rcs_ullage_time. // if any future active engine needs ullage, set ullage time.
    }
    print engine:name at (0,21+i).
    set i to i + 1.
}
print "Active engines: " + active_engines:length at (0,20).
set spool_time to spool_time / active_engines:length. // average the spoolup time of the active engines
print " Engine Spool Up Time: " + spool_time at (0,10).
local half_dv_time to calc_burn_time(total_dv/2).
print "         Half Dv Time: " + half_dv_time at (0,17).
lock t_burn to (nd:eta - half_dv_time).

// ############## Warp calculation ##############
until t_burn <= setup_time {
	local maxwarp to 5.
	if t_burn < 40000   { set maxwarp to 4. }
	if t_burn < 4000   { set maxwarp to 3. }
	if t_burn < 200 + setup_time    { set maxwarp to 2. }
	if t_burn < 50 + setup_time     { set maxwarp to 1. }
	if t_burn < setup_time + 2  { set maxwarp to 0. }
    set warp to maxwarp.
    print "          Warp factor: " + warp at (20,7).
    print_burn_time(t_burn).
    wait 0.1. // to avoid overloading the CPU
}


// ############## RCS attitude alignment ##############
RCS ON. // enable the RCS
lock steering to nd:deltav.

//now we need to wait until the delta-v vector and ship's facing are aligned
print_status("Waiting for attitude alignment...                  ").
until vang(nd:deltav, ship:facing:vector) < 0.25 {
    print_burn_time(t_burn).
    print_alignment().
}


// ############## RCS spin up ##############
if rcs_spin_time > 0 {
    print_status("Waiting for spin up time...                ").
    until t_burn <= spool_time/2 + ullage_time + rcs_spin_time {
        print_burn_time(t_burn).
        print_alignment().
    }
    print_status("Spinning up...             ").
    set ship:control:roll to 1. // turn on the RCS
    until t_burn <= spool_time/2 + ullage_time{
        print_burn_time(t_burn).
        print_alignment().
    }
} else {
    until t_burn <= spool_time/2 + ullage_time {
        print_burn_time(t_burn).
        print_alignment().
    }
    print_status("Waiting for fuel stability...             ").
}


// ############## RCS ullage ##############

when check_fuel_stability() < 1 then{
    set ship:control:fore to 1. // turn on the RCS
    print_status("Fuel unstable, doing ullage...  ").
}


until t_burn <= spool_time/2 - burn_delay {
    print_burn_time(t_burn).
    print_alignment().
    wait 0.
}


// ############## Stage separation ##############
set ship:control:roll to 0. // turn off the RCS
set ship:control:fore to 0. // turn off the RCS

runpath("0:/burn", stop_burn_tol, half_dv_time).