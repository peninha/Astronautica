from astronautica import Orbit, Body, Plotter, Frame, Maneuver, Trajectory
import numpy as np
import pandas as pd

earth = Body(name="earth")


e1 = 0
rp1 = earth.radius_from_altitude(250)
i1 = 12
Omega01 = 45
omega01 = 0
theta01 = 0
t0_clock1 = 0

orbit_from = Orbit.from_elements(earth, e=e1, rp=rp1, Omega0=Omega01, i=i1, omega0=omega01, theta0=theta01, t0_clock=t0_clock1)


e2 = 0.05
rp2 = earth.radius_from_altitude(400)
Omega02 = 45
i2 = 12
omega02 = 20
theta02 = 10
t0_clock2 = 0

orbit_to = Orbit.from_elements(earth, e=e2, rp=rp2, Omega0=Omega02, i=i2, omega0=omega02, theta0=theta02, t0_clock=t0_clock2)

omega_rel = Maneuver.relative_periapsis_argument(orbit_from, orbit_to)
print(f"omega_rel = {omega_rel} deg")

theta0_rel = Maneuver.relative_theta0(orbit_from, orbit_to)
print(f"theta0_rel = {theta0_rel} deg")

t1 = 3700

theta_rel_t_clock = Maneuver.relative_theta_at_t_clock(orbit_from, orbit_to, t1)
print(f"theta_rel_t_clock = {theta_rel_t_clock} deg")

traj1 = Trajectory(orbit_from)
traj1.add_trajectory_position(0, t_clock=t1)
traj2 = Trajectory(orbit_to)
traj2.add_trajectory_position(0, t_clock=t1)



plotter = Plotter(plot3d=True)
plotter.plot_trajectories([traj1, traj2], frame="bodycentric")
