runpath("snd_stage:/roll.ks").
stage.
print "RCS Dv: " + rcs_dv at (3,11).
set t0 to time:seconds. lock throttle to 1. set b_dv to rcs_dv.
until b_dv >= tg_dv - 1.5 {set t1 to time:seconds. set b_dv to b_dv + (ship:thrust/ship:mass)*(t1-t0). set t0 to t1. print "Target Dv: " + tg_dv at (0,12). print "Burned Dv: " + b_dv at (0,13).}
lock throttle to 0.
print "Vai FILHÃÃÃO!".
TOGGLE AG5.