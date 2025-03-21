CLEARSCREEN.
// Using Delta-V to calculate the burn

// ############## Parameters ##############
set nd to nextnode.
set setup_time to 90. // setup time for the burn
set rcs_ullage_time to 3. // ullage time for the burn
set rcs_isp to 136.5. // ISP for the RCS burn
set g0 to constant():G0. // gravity constant
set stop_burn_tol to 0.1. // tolerance for the burn (delta_v before reaching target to stop the burn)


// ############## Functions ##############
function print_status {
    parameter info.
    print "status: " +info at (5, 5).
}

function print_burn_time {
    print " Burn time: " + ROUND(-1*nd:eta,1) + " s      " at (0,7).
}

function print_alignment {
    set align to vang(nd:deltav, ship:facing:vector).
    if align < 0.25 print "         Alignment OK: " + ROUND(align,2) + "deg    " AT(0,9).
    else print "            Alignment: " + ROUND(align,2) + "deg    " AT(0,9).
}

global active_engines is list().
function check_fuel_stability{
    set ready to 1.
    for engine in active_engines {
        set ready to ready * engine:fuelstability.
    }
    return ready.
}

print "######### Next Node Burning Program #########" at(0, 3).
print_status("Warping to node...                        ").


// ############## Active engines calculations ##############

set ullage_time to 0.
set spool_time to 0.
list engines in MyEngines.
set active_engines to list().
for engine in MyEngines {
	if engine:ignition  {
        // add the spoolup time of the active engines
        set spool_time to spool_time + 0.414213 * ln(max(1.1,sqrt(engine:drymass*engine:maxthrust^2))).
        active_engines:add(engine). // add the engine to the list of active engines
        if engine:ullage {
            set ullage_time to rcs_ullage_time. // if any active engine needs ullage, set ullage time.
        }
	}
}
set spool_time to spool_time / active_engines:length. // average the spoolup time of the active engines
print " Engine Spool Up Time: " + spool_time at (0,10).


// ############## Warp calculation ##############
until nd:eta <= setup_time {
	set maxwarp to 5.
    set t_burn to nd:eta.
	if t_burn < 50000   { set maxwarp to 4. }
	if t_burn < 5000   { set maxwarp to 3. }
	if t_burn < 500    { set maxwarp to 2. }
	if t_burn < 200     { set maxwarp to 1. }
	if t_burn < setup_time + 5   { set maxwarp to 0. }
    set warp to maxwarp.
    print "          Warp factor: " + warp at (20,7).
    print_burn_time().
    wait 1. // to avoid overloading the CPU
}


// ############## RCS attitude alignment ##############
RCS ON. // enable the RCS
lock steering to nd:deltav.

//now we need to wait until the delta-v vector and ship's facing are aligned
print_status("Waiting for attitude alignment...                  ").
until vang(nd:deltav, ship:facing:vector) < 0.25 {
    print_burn_time().
    print_alignment().
}
print_status("Waiting for ullage time...                ").

until nd:eta <= spool_time + ullage_time {
    print_burn_time().
    print_alignment().
}
print_status("Waiting for fuel stability...             ").


// ############## RCS ullage ##############
set rcs_dv to 0.
if check_fuel_stability() < 1 {
    set m0 to ship:mass.
    set ship:control:fore to 1. // turn on the RCS
    set ve to rcs_isp * g0. // effective exhaust velocity for RCS ullage
    until false{
        if check_fuel_stability() = 1{
            print_status("Fuel stable, waiting for burning time...  ").
            if nd:eta <= spool_time/2{ // half the spool time to average the thrust before and after the burn.
                break.
            }
        }
        set mf to ship:mass.
        set rcs_dv to rcs_dv + ve * ln(m0/mf).
        set m0 to mf.
        print_burn_time().
        print_alignment().
        print "            RCS Dv: " + rcs_dv at (0,11).

        wait 0.
    }
    set ship:control:fore to 0. // turn off the RCS
}


// ############## Engines spooling up ##############
set t0 to time:seconds.
set ship:control:pilotmainthrottle to 1. // turn on the engines
print_status("Engines spooling up...                    ").


// ############## Engines burn ##############
when nd:eta <= -spool_time/2 then {
    print_status("Engines burning...                        ").
}

set burn_dv to rcs_dv. // include RCS ullage delta-v for the burn
set target_dv to nd:deltav:mag.

until burn_dv >= target_dv - stop_burn_tol {
    local t1 to time:seconds.  
    set burn_dv to burn_dv + (ship:thrust/ship:mass)*(t1-t0).
    set t0 to t1.

    print_burn_time().
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



