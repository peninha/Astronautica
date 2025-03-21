clearscreen. print "running...".
wait until ship:mass < 0.8. set t0 to time:seconds. lock throttle to 1. set b_dv to 0. set tg_dv to nextnode:deltav:mag.
until b_dv >= tg_dv - 1.7 {
set t1 to time:seconds. set b_dv to b_dv + (ship:thrust/ship:mass)*(t1-t0). set t0 to t1. print "Target Dv: " + tg_dv at (0,12). print "Burned Dv: " + b_dv at (0,13).
}
lock throttle to 0. print "Finished!".