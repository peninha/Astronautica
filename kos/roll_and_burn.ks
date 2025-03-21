CLEARSCREEN.
// Using Delta-V to calculate the burn

// ############## Parameters ##############
local nd to nextnode.
local setup_time to 120. // setup time for the burn
local rcs_spin_time to 30. // time to spin the ship before the burn
local rcs_ullage_time to 2. // ullage time for the burn
//local rcs_isp to 136.5. // ISP for the RCS burn
local g0 to constant():G0. // gravity constant
local stop_burn_tol to 0.4. // tolerance for the burn (delta_v before reaching target to stop the burn)
//local total_dv to 3224. // total delta-v for the burn
//local stage_mass to 0.420. // mass of the spin burn stage
local stage_mass to 0.743. // mass of the spin burn stage
local burn_delay to 2.2. // delay for the burn
local target_dv to nd:deltav:mag.

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
local half_dv_time to calc_burn_time(target_dv/2).
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
    //set ship:control:roll to 0. // turn off the RCS
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

local rcs_dv to 0.
until t_burn <= spool_time/2 - burn_delay {
    //local m0 to ship:mass.
    //local ve to rcs_isp * g0. // effective exhaust velocity for RCS ullage
    //local mf to ship:mass.
    //set rcs_dv to rcs_dv + ve * ln(m0/mf).
    //set m0 to mf.

    print_burn_time(t_burn).
    print_alignment().
    //print "            RCS Dv: " + rcs_dv at (0,11).
    wait 0.
}

// ############## Stage separation ##############
set ship:control:roll to 0. // turn off the RCS
set ship:control:fore to 0. // turn off the RCS
stage.

// ############## Engines spooling up ##############
local t0 to time:seconds.
print "    burn initial mass: " + ship:mass at (0,15).
print "             dry mass: " + ship:drymass at (0,16).

set ship:control:pilotmainthrottle to 1. // turn on the engines
print_status("Engines spooling up...                    ").


// ############## Engines burn ##############
when t_burn <= (-spool_time/2) - burn_delay then {
    print_status("Engines burning...                        ").
}

local burn_dv to rcs_dv. // include RCS ullage delta-v for the burn


until burn_dv >= target_dv - stop_burn_tol {
    local t1 to time:seconds.  
    set burn_dv to burn_dv + (ship:thrust/ship:mass)*(t1-t0).
    set t0 to t1.

    print_burn_time(t_burn).
    print_alignment().
    print "            Target Dv: " + target_dv at (0,12).
    print "      Total Burned Dv: " + burn_dv at (0,13).
    print "           deltav:mag: " + nd:deltav:mag at (0,14).

    wait 0.
}

set ship:control:pilotmainthrottle to 0. // turn off the engines
print_status("Finished burn successfully!                ").
unlock all.

//we no longer need the maneuver node
//remove nd.