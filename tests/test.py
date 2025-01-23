from astronautica import Orbit, Body, Plotter, Frame, Maneuver, Trajectory
import numpy as np
import pandas as pd

earth = Body(name="spherical_earth")


e1 = 0
rp1 = earth.radius_from_altitude(250)
i1 = 12
Omega01 = 45
omega01 = 0
theta01 = 0
t0_clock1 = 0

orbit_from = Orbit.from_elements(earth, e=e1, rp=rp1, Omega0=Omega01, i=i1, omega0=omega01, theta0=theta01, t0_clock=t0_clock1)

rp2 = earth.radius_from_altitude(1400)
ra2 = earth.radius_from_altitude(3000)
Omega02 = 45
i2 = 12
omega02 = 20
theta02 = 10
t0_clock2 = 0

orbit_to = Orbit.from_elements(earth, rp=rp2, ra=ra2, Omega0=Omega02, i=i2, omega0=omega02, theta0=theta02, t0_clock=t0_clock2)

loop_period = Maneuver.relative_theta_loop_period(orbit_from, orbit_to, from_t_clock=2000)
print(f"loop_period = {loop_period} s")

theta_rel = -120
rel_t_clock = Maneuver.find_t_clock_for_relative_theta(orbit_from, orbit_to, theta_rel, from_t_clock=2000, max_iterations=1000, tol=1e-6)
print(f"rel_t_clock = {rel_t_clock} s")

theta_rel = Maneuver.relative_theta_at_t_clock(orbit_from, orbit_to, rel_t_clock)
print(f"theta_rel = {theta_rel} deg")

traj1 = Trajectory(orbit_from)
traj1.add_trajectory_position(0, t_clock=rel_t_clock)
traj2 = Trajectory(orbit_to)
traj2.add_trajectory_position(0, t_clock=rel_t_clock)


plotter = Plotter(plot3d=False)
plotter.plot_trajectories([traj1, traj2], frame="bodycentric",
                           time_step=60,
                           orbits=True,
                           orbit_arrows=False,
                           points=True,
                           positions=True,
                           velocities=True)
