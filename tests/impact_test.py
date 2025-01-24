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
i1 = moon_orbit.i - 0
Omega01 = moon_orbit.Omega0 + 0
omega01 = 0
theta01 = 0
t0_clock1 = 0
parking_orbit = Orbit.from_elements(earth, e=e1, rp=rp1, Omega0=Omega01, i=i1, omega0=omega01, theta0=theta01, t0_clock=t0_clock1)

theta_rel = 60
from_t_clock = 800
t_clock_burn = Maneuver.find_t_clock_for_relative_theta(parking_orbit, moon_orbit, theta_rel, from_t_clock=from_t_clock, max_iter=1000, tol=1e-9)
print(f"t_clock_burn = {t_clock_burn} s")

########################
# delta_v_target = 1.5 #
########################
delta_t = 3900
t_clock_impact = t_clock_burn + delta_t

r_vec_from, v_vec_from = parking_orbit.state_vectors_at_t_clock(t_clock_burn)
r_vec_impact, v_vec_impact = moon_orbit.state_vectors_at_t_clock(t_clock_impact)

impact_orbit = Orbit.from_2_vectors_and_delta_time(earth, r_vec_from, r_vec_impact, delta_t, t0_clock=t_clock_impact)

theta_burn = parking_orbit.theta_at_t_clock(t_clock_burn)
print(f"theta_burn = {theta_burn} deg")
burn, _ , delta_v = Maneuver.orbit_change_maneuver_at_theta(parking_orbit, impact_orbit, theta_burn)
print(f"delta_v = {delta_v} m/s")


traj1 = Trajectory(parking_orbit, orbit_name="Parking orbit", position_name="Probe initial position")
traj1.add_maneuver(0, burn, name="burn", position_name="Probe burn position")
traj1.add_trajectory_position(1, t_clock=t_clock_impact, name="Probe impact position")

traj2 = Trajectory(moon_orbit, orbit_name="Moon orbit", position_name="Moon initial position")
traj2.add_trajectory_position(0, t_clock=t_clock_burn, name="Moon at burn")
traj2.add_trajectory_position(0, t_clock=t_clock_impact, name="Moon at impact")

plotter = Plotter(plot3d=True)
plotter.plot_trajectories([traj1, traj2], frame="perifocal",
                           time_step=60,
                           orbits=False,
                           orbit_arrows=True)


