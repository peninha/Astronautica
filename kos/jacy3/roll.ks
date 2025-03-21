clearscreen. print "partiu...".
set set_t to 100, spin_t to 40, u_t to 7, dl to 0, dv_adj to -2, rcs_thr to 0.08, spool to 0.16.

set nd to nextnode, tg_dv to nd:deltav:mag+dv_adj.
set h_dv_t to (0.519 * 231 * constant():G0 / 21.28) * (1 - constant():e^(-tg_dv/2 / (231*constant():G0))).
lock t_burn to (nd:eta - h_dv_t).
set e to ship:engines[1]. e:activate().

wait until t_burn <= set_t.
RCS ON. lock steering to nd:deltav.

wait until t_burn <= spool/2 + u_t + spin_t.
TOGGLE AG6. set ship:control:roll to 1.

wait until t_burn <= spool/2 + u_t.
set t0 to time:seconds.
when e:fuelstability < 1 then{set ship:control:roll to 0. unlock steering. set ship:control:fore to 1. set t0 to time:seconds.}

set rcs_dv to 0.
until t_burn <= spool/2 - dl {
    if ship:control:fore = 1{
        local t1 to time:seconds.  
        set rcs_dv to rcs_dv + (rcs_thr/ship:mass)*(t1-t0).
        set t0 to t1.
    }
}
wait until t_burn <= spool/2 - dl.
set ship:control:roll to 0. set ship:control:fore to 0.