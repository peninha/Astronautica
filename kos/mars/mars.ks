CLEARSCREEN.
// Using Delta-V to calculate the burn

// ############## Parameters ##############
local nd to nextnode.
local setup_time to 100. // setup time for the burn
local rcs_spin_time to 50. // time to spin the ship before the burn
local rcs_ullage_time to 10. // ullage time for the burn
//local rcs_isp to 144.4. // ISP for the RCS burn
//local rcs_thr to 0.09226. // thrust for the RCS burn
local g0 to constant():G0. // gravity constant
local stop_burn_tol to 1.5. // tolerance for the burn (delta_v before reaching target to stop the burn)
//local stage_mass to 1.043. // mass of the spin burn stage
local burn_delay to 0. // delay for the burn
local dv_adjust to 0. // delta-v adjustment for the burn
local target_dv to nd:deltav:mag + dv_adjust.
local stage2_dv to 1183.
local half_dv_time to 47.
local stage1_dv to target_dv - stage2_dv - 5.5.

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
    if align < 0.25 print " Alignment OK: " + ROUND(align,2) + "deg    " AT(0,8).
    else print " Alignment: " + ROUND(align,2) + "deg    " AT(0,8).
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

print "######### Roll and Burn Program #########" at(0, 3).
print_status("Warping to node...                        ").

// ############## Active engines calculations ##############

local ullage_time to 0.
local spool_time to 0.
list engines in engine_list.
local i to 0.
for engine in engine_list {
    if engine:ignition{ // Add active engines to the list.
        active_engines:add(engine). // add the engine to the list
	    set spool_time to spool_time + calc_spool_time(engine).
        if engine:ullage {
            set ullage_time to rcs_ullage_time. // if any future active engine needs ullage, set ullage time.
        }
        print engine:name at (0,30+i).
        set i to i + 1.
    }
}
print "Active engines: " + active_engines:length at (0,29).
set spool_time to spool_time / active_engines:length. // average the spoolup time of the active engines
print " Engine Spool Up Time: " + ROUND(spool_time,2) at (0,10).
//local half_dv_time to calc_burn_time(target_dv/2).
print "         Half Dv Time: " + ROUND(half_dv_time,2) at (0,11).
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
    set ship:control:roll to 1. // turn on roll
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
local ullage_motors is list().
local rcs_thrust to 0.
local rcs_isp to 0.
LIST RCS IN RCS_list.
for rcs_motor in RCS_list {
    if rcs_motor:enabled and rcs_motor:foreenabled{ // check only enabled foreward facing nozzles
        ullage_motors:add(rcs_motor). // add the engine to the list
	    set rcs_thrust to rcs_thrust + rcs_motor:availablethrust.
        set rcs_isp to rcs_isp + rcs_motor:visp * rcs_motor:availablethrust. // isp weighted by thrust
    }
}
set rcs_isp to rcs_isp / rcs_thrust. // average isp by thrust
print "            RCS ISP: " + ROUND(rcs_isp,3) at (0,12).
print "            RCS Thrust: " + ROUND(rcs_thrust,5) at (0,13).

local t0 to time:seconds.
when check_fuel_stability() < 1 then{
    set ship:control:roll to 0. // turn off roll
    unlock steering.
    set ship:control:fore to 1. // turn on the ullage
    set t0 to time:seconds.
    print_status("Fuel unstable, doing ullage...  ").
}

local rcs_dv to 0.
until t_burn <= spool_time/2 - burn_delay {
    if ship:control:fore = 1{
        local t1 to time:seconds.  
        set rcs_dv to rcs_dv + (rcs_thrust/ship:mass)*(t1-t0).
        set t0 to t1.
    }
    print_burn_time(t_burn).
    print_alignment().
    print "            RCS Dv: " + ROUND(rcs_dv,2) at (0,14).
    wait 0.
}

// ############## Stage separation ##############
set ship:control:roll to 0. // turn off the RCS
set ship:control:fore to 0. // turn off the RCS
stage.

// ############## Engines spooling up ##############
set t0 to time:seconds.
print "    burn initial mass: " + ROUND(ship:mass,6) at (0,15).

set ship:control:pilotmainthrottle to 1. // turn on the engines
print_status("Engines spooling up...                    ").


// ############## Engines burn ##############
when t_burn <= (-spool_time/2) - burn_delay then {
    print_status("Engines burning...                        ").
}

local burn_dv to rcs_dv. // include RCS ullage delta-v for the burn

until burn_dv >= stage1_dv - stop_burn_tol {
    local t1 to time:seconds.  
    set burn_dv to burn_dv + (ship:thrust/ship:mass)*(t1-t0).
    set t0 to t1.

    print_burn_time(t_burn).
    print_alignment().
    print "    STAGE 1" at (0,17).
    print "            Target Dv: " + ROUND(stage1_dv,2) at (0,18).
    print "            Burned Dv: " + ROUND(burn_dv,2) at (0,19).
    print "      Total Target Dv: " + ROUND(target_dv,2) at (0,25).
    print "      Total Burned Dv: " + ROUND(burn_dv,2) at (0,26).
    wait 0.
}
set ship:control:pilotmainthrottle to 0. // turn off the engines
stage.
set ship:control:pilotmainthrottle to 1. // turn on the engines
set t0 to time:seconds.
local burn_dv_prev to burn_dv.
local end to false.
until end {
    local t1 to time:seconds.
    set burn_dv to burn_dv + (ship:thrust/ship:mass)*(t1-t0).
    set t0 to t1.

    print_burn_time(t_burn).
    print_alignment().
    print "    STAGE 2" at (0,21).
    print "            Target Dv: " + ROUND(stage2_dv,2) at (0,22).
    print "            Burned Dv: " + ROUND(burn_dv - burn_dv_prev,2) at (0,23).
    print "      Total Target Dv: " + ROUND(target_dv,2) at (0,25).
    print "      Total Burned Dv: " + ROUND(burn_dv,2) at (0,26).
    wait 0.
    if ship:thrust = 0 {
        set end to true.
    }
}

set ship:control:pilotmainthrottle to 0. // turn off the engines
print_status("Finished burn successfully!                ").
unlock all.

//we no longer need the maneuver node
//remove nd.