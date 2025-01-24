from astronautica import Orbit, Body, Plotter, Frame, Maneuver, Trajectory
import numpy as np
import pandas as pd

earth = Body(name="spherical_earth")

rp = earth.radius_from_altitude(5000)
ra = earth.radius_from_altitude(8000)
i = 18
Omega0 = 25
omega0 = -45
theta0 = 52
t0_clock = 0
moon_orbit = Orbit.from_elements(earth, rp=rp, ra=ra, Omega0=Omega0, i=i, omega0=omega0, theta0=theta0, t0_clock=t0_clock)

e1 = 0
rp1 = earth.radius_from_altitude(250)
i1 = moon_orbit.i - 10
Omega01 = moon_orbit.Omega0 + 10
omega01 = 0
theta01 = 0
t0_clock1 = 0
parking_orbit = Orbit.from_elements(earth, e=e1, rp=rp1, Omega0=Omega01, i=i1, omega0=omega01, theta0=theta01, t0_clock=t0_clock1)

theta_rel = 60
from_t_clock = 800
delta_v_target = 2.5
delta_t_min = 3000
delta_t_max = 4000
delta_t_guess = 3500

result = Maneuver.impact_maneuver(parking_orbit, moon_orbit, theta_rel, delta_v_target, delta_t_guess, delta_t_min, delta_t_max, from_t_clock=from_t_clock, tol=1e-6, max_iter=1000)

burn = result["maneuver"]
delta_t = result["delta_t"]
delta_v = result["delta_v"]
t_clock_burn = result["t_clock_burn"]
t_clock_impact = t_clock_burn + delta_t

print("delta_t: ", delta_t)
print("delta_v: ", delta_v)
print("t_clock_burn: ", t_clock_burn)
print("t_clock_impact: ", t_clock_impact)

traj1 = Trajectory(parking_orbit, orbit_name="Parking orbit", position_name="Probe initial position")
traj1.add_maneuver(0, burn, name="Translunar injection")
traj1.add_trajectory_position(1, t_clock=t_clock_impact, name="Probe impact position")

traj2 = Trajectory(moon_orbit, orbit_name="Moon orbit", position_name="Moon initial position")
traj2.add_trajectory_position(0, t_clock=t_clock_burn, name="Moon at burn")
traj2.add_trajectory_position(0, t_clock=t_clock_impact, name="Moon at impact")

plotter = Plotter(plot3d=True)
plotter.plot_trajectories([traj1, traj2], frame="perifocal",
                           time_step=60,
                           orbits=False)

