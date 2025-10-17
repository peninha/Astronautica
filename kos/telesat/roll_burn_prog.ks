CLEARSCREEN.
// Using Delta-V to calculate the burn

// ############## Parameters ##############
parameter angular_vel to 60.00. // intended angular velocity in degrees per second
parameter align_tol to 0.10. // alignment tolerance in degrees
local sample_time to 30.0. // time to collect samples for wobble calculation
local sample_interval to 0.1. // time between samples for wobble calculation

local spool_time to 0.
local rcs_dv to 0.
local rcs_thrust to 0.
local rcs_isp to 0.
lock t_burn to 0.

// ############## Functions ##############
function print_status {
    parameter info.
    print "status: " +info at (5, 5).
}

function print_burn_time {
    parameter t_burn.
    print "            Burn time: " + ROUND(t_burn,1) + " s      " at (0,16).
}

function print_alignment {
    parameter prograde_vec to ship:velocity:orbit.
    local align to vang(prograde_vec, ship:facing:vector).
    if align < align_tol print "         Alignment OK: " + ROUND(align,2) + " deg    " at(0,7).
    else print "            Alignment: " + ROUND(align,2) + " deg    " at(0,7).
}

function print_angular_velocity {
    print "     Angular velocity: " + ROUND(ship:angularvel:mag * constant:radtodeg, 2) + " deg/s    " at (0,9). 
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

function calc_wobble{
    parameter sample_time is 30.0.
    parameter sample_interval is 0.1.

    local samples is round(sample_time / sample_interval).
    local dt is sample_interval.

    // Lista para armazenar as amostras de angularvel
    local facing_vecs is list().

    local avg is V(0,0,0).
    local i is 0.
    until i = samples {
        WAIT dt.
        facing_vecs:add(ship:facing:vector).
        set avg to avg + facing_vecs[i].
        set i to i + 1.
    }
    set avg to avg / samples.

    // Calculate wobbleness (average angle between each facing vector and the average facing vector)
    set total_angle to 0.
    for sample IN facing_vecs {
        set angle_diff to vang(avg, sample).
        set total_angle to total_angle + angle_diff.
    }
    set average_angle to total_angle / samples.
    return average_angle.
}

print "######### Roll and Burn Prograde #########" at(0, 3).

print "  Target angular vel.: " + ROUND(angular_vel, 2) + " deg/s    " at(0,8).

// ############## Active engines calculations ##############
list engines in engine_list.
local i to 0.
for engine in engine_list {
    if engine:ignition{ // Add active engines to the list.
        active_engines:add(engine). // add the engine to the list
	    set spool_time to spool_time + calc_spool_time(engine).
        print engine:name at (0,20+i).
        set i to i + 1.
    }
}
print "Active engines: " + active_engines:length at (0,19).
set spool_time to spool_time / active_engines:length. // average the spoolup time of the active engines
print " Engine Spool Up Time: " + ROUND(spool_time, 2) + " s    " at (0,11).

// ############## RCS attitude alignment ##############
RCS ON. // enable the RCS
lock steering to ship:velocity:orbit.

//now we need to wait until the delta-v vector and ship's facing are aligned
print_status("Waiting for attitude alignment...       ").
until vang(ship:velocity:orbit, ship:facing:vector) < align_tol {
    print_burn_time(t_burn).
    print_alignment().
    print_angular_velocity().
}

// ############## RCS spin up ##############
until ship:angularvel:mag * constant:radtodeg >= angular_vel {
    print_status("Spinning up...                          ").
    set ship:control:roll to 1. // turn on roll
    print_burn_time(t_burn).
    print_alignment().
    print_angular_velocity().
}

// ############## RCS ullage ##############
local ullage_motors is list().
LIST RCS IN RCS_list.
for rcs_motor in RCS_list {
    if rcs_motor:enabled and rcs_motor:foreenabled{ // check only enabled foreward facing nozzles
        ullage_motors:add(rcs_motor). // add the engine to the list
	    set rcs_thrust to rcs_thrust + rcs_motor:availablethrust.
        set rcs_isp to rcs_isp + rcs_motor:visp * rcs_motor:availablethrust. // isp weighted by thrust
    }
}
set rcs_isp to rcs_isp / rcs_thrust. // average isp by thrust
print "              RCS ISP: " + ROUND(rcs_isp, 2) + " s    " at (0,12).
print "           RCS Thrust: " + ROUND(rcs_thrust, 2) + " kN    " at (0,13).
print "               RCS Dv: " + ROUND(rcs_dv,2) + " m/s    " at (0,14).

local prograde_vec to ship:velocity:orbit.
unlock steering.
set ship:control:roll to 0. // turn off roll

local t0 to time:seconds.
if check_fuel_stability() < 1 {
    print_status("Fuel unstable, doing ullage...          ").
    set ship:control:fore to 1. // turn on the ullage
    set t0 to time:seconds.

    until check_fuel_stability() = 1 {
        local t1 to time:seconds.  
        set rcs_dv to rcs_dv + (rcs_thrust/ship:mass)*(t1-t0).
        set t0 to t1.
        print_burn_time(t_burn).
        print_alignment(prograde_vec).
        print_angular_velocity().
        print "               RCS Dv: " + ROUND(rcs_dv,2) + " m/s    " at (0,14).
        wait 0.
    }
}

// ############## Stage separation ##############
set ship:control:fore to 0. // turn off the RCS
stage.

// ############## Engines spooling up ##############
set t0_burn to time:seconds.
lock t_burn to time:seconds - t0_burn.
set ship:control:pilotmainthrottle to 1. // turn on the engines
set t0 to time:seconds.
print_status("Engines spooling up...                  ").

// ############## Engines burn ##############
when t_burn >= (spool_time/2) then {
    print_status("Engines burning...                      ").
}
local burn_dv to rcs_dv. // include RCS ullage delta-v for the burn

local end to false.
until end {
    local t1 to time:seconds.  
    set burn_dv to burn_dv + (ship:thrust/ship:mass)*(t1-t0).
    set t0 to t1.
    print_burn_time(t_burn).
    print_alignment(prograde_vec).
    print_angular_velocity().
    print "            Burned Dv: " + ROUND(burn_dv,2) + " m/s    " at (0,17).
    wait 0.
    if ship:thrust = 0 and t_burn > (spool_time + 1) {
        set end to true.
    }
}

set ship:control:pilotmainthrottle to 0. // turn off the engines
print_status("Finished burn, waiting to stage...      ").
print "                                " at (0,19).
print "                                " at (0,20).
print "                                " at (0,21).
print "                                " at (0,22).
print "                                " at (0,23).
print "                                " at (0,24).
print "                                " at (0,25).
wait 2.
stage.
print_status("Staged. Checking final stability...     ").
wait 1.

print "Getting samples... wait " + sample_time + " seconds." at (0,19).
local wobble is calc_wobble(sample_time, sample_interval).
print "wobble: " + wobble + " deg" at (0,20).
print_status("Finished program successfully!          ").
unlock all.