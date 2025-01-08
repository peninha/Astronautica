from astronautica import Orbit, Body, Plotter, Frame, Maneuver, Trajectory
import numpy as np

earth = Body(name="spherical_earth")

rp = 12000
ra = 35000
theta0 = 90
theta_burn = 150
t0_clock = 0
Omega0 = 40
i = 40
omega0 = 0

orbit_from = Orbit.from_elements(earth, rp=rp, ra=ra, Omega0=Omega0, i=i, omega0=omega0, theta0=theta0, t0_clock=t0_clock)
traj1 = Trajectory(orbit_from, orbit_name="Starting orbit", position_name="Initial position")
traj1.add_trajectory_position(0, theta=theta_burn, name="Transfer position")

rp = 8000
ra = 10000
omega0 = 30
theta0 = 0

orbit_to = Orbit.from_elements(earth, rp=rp, ra=ra, Omega0=Omega0, i=i, omega0=omega0, theta0=theta0, t0_clock=t0_clock, name="Final orbit", position_name="Arrival position")

maneuver1, maneuver2, delta_v = Maneuver.hohmann_transfer(orbit_from, orbit_to, theta_burn)

print("delta_v: ", delta_v)

traj1 = Trajectory(orbit_from, orbit_name="Starting orbit", position_name="Starting position")
traj1.add_maneuver(0, maneuver1)
traj1.add_maneuver(1, maneuver2)
traj1.add_trajectory_position(2, t_clock=traj1.get_trajectory_position(2, "last")["t_clock"] + traj1.orbits[2]["orbit"].T/2, name="End position")

plotter = Plotter(plot3d=True)
plotter.plot_trajectory(traj1, frame="bodycentric", orbits=True, velocities=True, positions=True, maneuvers=True, groundtrack=False, v_scale=3)
