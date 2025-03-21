CLEARSCREEN.
// Using Delta-V to calculate the burn

// ############## Parameters ##############
set nd to nextnode.
set setup_time to 25. // setup time for the burn
set rcs_time to 3. // ullage time for the burn
set rcs_isp to 187.2. // ISP for the RCS burn
set g0 to constant():G0. // gravity constant
set stop_burn_tol to 1. // tolerance for the burn (delta_v before reaching target to stop the burn)
set delta_v_vec to nd:deltav. // burn vector


// ############## Functions ##############
function print_status {
    parameter info.
    print "status: " +info at (5, 5).
}

function print_burn_time {
    print " Burn time: " + ROUND(-1*nd:eta,1) + " s      " at (0,7).
}

function print_alignment {
    set align to vang(delta_v_vec, ship:facing:vector).
    if align < 0.25 print "         Alignment OK: " + ROUND(align,2) + "deg    " AT(0,9).
    else print "            Alignment: " + ROUND(align,2) + "deg    " AT(0,9).
}

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
set ullage to 0.
set spooltime to 0.
set engines_isp to 0.
list engines in MyEngines.
set active_engines to list().
for engine in MyEngines {
	if engine:ignition  {
        // use the higher spoolup time of the active engines
        set spooltime to max(spooltime, 0.414213 * ln(max(1.1,sqrt(engine:drymass*engine:maxthrust^2)))).
        set engines_isp to engines_isp + engine:visp. // add the ISP of the active engines
        active_engines:add(engine). // add the engine to the list of active engines
        if engine:ullage {
            set ullage to rcs_time. // if any active engine needs ullage, set ullage time.
        }
	}
    set engines_isp to engines_isp / active_engines:length. // average the ISP of the active engines
}
print " Engine Spool Up Time: " + spooltime at (0,10).


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
RCS ON.
lock steering to delta_v_vec.

//now we need to wait until the burn vector and ship's facing are aligned
print_status("Waiting for attitude alignment...                  ").
until vang(delta_v_vec, ship:facing:vector) < 0.25 {
    print_burn_time().
    print_alignment().
}
print_status("Waiting for ullage time...                ").


// ############## Wait for ullage time ##############
until nd:eta <= spooltime + ullage {
    print_burn_time().
    print_alignment().
}
print_status("Waiting for fuel stability...             ").


// ############## RCS ullage ##############
set m0 to ship:mass. // initial mass before ullage
set ship:control:fore to 1. // turn on the RCS
until check_fuel_stability() = 1 {
    print_burn_time().
    print_alignment().
}
print_status("Fuel stable, waiting for burning time...  ").

until nd:eta <= spooltime {
    print_burn_time().
    print_alignment().
}
print_status("Engines spooling up...                    ").
set ship:control:fore to 0. // turn off the RCS
set ve to rcs_isp * g0. // effective exhaust velocity for RCS ullage
set mf to ship:mass. // final mass after RCS ullage
set rcs_dv to ve * ln(m0/mf). // delta-v for RCS ullage
print "          RCS delta-v: " + rcs_dv at (0,11).


// ############## Engines burn ##############
set ship_v0_vec to ship:velocity:orbit. // initial velocity vector
set m0 to ship:mass. // initial mass before burn
print "         Initial Mass: " + m0 at (0,15).
set ship:control:pilotmainthrottle to 1. // turn on the engines
when nd:eta <= 0 then {
    print_status("Engines burning...                        ").
}

//set burn_dv to rcs_dv. // include RCS ullage delta-v for the burn
//set avg_isp to 0.
//set numsamples to 0.
//lock mf to ship:mass. // track the final mass

//set target_dv to nd:deltav:mag.
//set target_mass to m0 / (constant:e ^ ((target_dv - rcs_dv) / (engines_isp * g0))).


local maneuver_dv to nextnode:deltav:mag.
local spent_dv to 0.
local t0 to time:seconds.

// ignition

until false {
  local t1 to time:seconds.  
  local dt to t1-t0.
  local acc to (ship:thrust/ship:mass)*dt.

  set spent_dv to spent_dv + acc. 
  set t0 to t1.

  if spent_dv >= maneuver_dv { break. }

  wait 0.
}



set ship:control:pilotmainthrottle to 0. // turn off the engines
print_status("Finished burn successfully!                ").
unlock all.

//we no longer need the maneuver node
//remove nd.



