parameter stop_burn_tol, half_dv_time.

set nd to nextnode.

lock t_burn to (nd:eta - half_dv_time).

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

stage.

// ############## Engines spooling up ##############
local t0 to time:seconds.
print "    burn initial mass: " + ship:mass at (0,15).
print "             dry mass: " + ship:drymass at (0,16).

set ship:control:pilotmainthrottle to 1. // turn on the engines
print_status("Engines spooling up...                    ").


local burn_dv to 0.
local target_dv to nd:deltav:mag.

until burn_dv >= target_dv - stop_burn_tol {
    local t1 to time:seconds.  
    set burn_dv to burn_dv + (ship:thrust/ship:mass)*(t1-t0).
    set t0 to t1.

    print_burn_time(t_burn).
    print_alignment().
    print "            Target Dv: " + target_dv at (0,12).
    print "      Total Burned Dv: " + burn_dv at (0,13).

    wait 0.
}

set ship:control:pilotmainthrottle to 0.
print_status("Finished burn successfully!                ").
unlock all.
