parameter set_t to 120, spin_t to 30, u_t to 2, delay to 2.2.

clearscreen.
print "running...".

local spool to 0.16.
ship:engines[0]:activate().

local h_dv_t to (0.42 * 231 * constant():G0 / 21.28) * (1 - constant():e^(-nextnode:deltav:mag/2 / (231*constant():G0))).
lock t_burn to (nextnode:eta - h_dv_t).

until t_burn <= set_t {
	local maxwarp to 5.
	if t_burn < 40000   { set maxwarp to 4. }
	if t_burn < 4000   { set maxwarp to 3. }
	if t_burn < 200 + set_t    { set maxwarp to 2. }
	if t_burn < 50 + set_t     { set maxwarp to 1. }
	if t_burn < set_t + 2  { set maxwarp to 0. }
    set warp to maxwarp.
}

RCS ON.
lock steering to nextnode:deltav.

wait until t_burn <= spool/2 + u_t + spin_t.
set ship:control:roll to 1.


wait until t_burn <= spool/2 + u_t.
when ship:engines[0]:fuelstability < 1 then{
    set ship:control:fore to 1.
}

wait until t_burn <= spool/2 - delay.
set ship:control:roll to 0.
set ship:control:fore to 0.
stage.