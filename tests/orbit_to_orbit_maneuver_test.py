from astronautica import Orbit, Body, Plotter, Frame, Maneuver, Trajectory
import numpy as np

earth = Body(name="spherical_earth")

rp = 9000
ra = 35000
theta0 = -180
theta_burn = -150
t0_clock = 0
Omega0 = 0
i = 60
omega0 = -40

orbit_from = Orbit.from_elements(earth, rp=rp, ra=ra, Omega0=Omega0, i=i, omega0=omega0, theta0=theta0, t0_clock=t0_clock)

rp = 8000
ra = 8000
theta_arrival = 0

orbit_to = Orbit.from_elements(earth, rp=rp, ra=ra, Omega0=Omega0, i=i, omega0=omega0, theta0=theta_arrival, t0_clock=t0_clock, name="Final orbit", position_name="Arrival position")

maneuver1, maneuver2, delta_t, delta_v = Maneuver.orbit_pos_to_orbit_pos_maneuver(orbit_from, orbit_to, theta_burn, theta_arrival)

traj1 = Trajectory(orbit_from, orbit_name="Starting orbit", position_name="Starting position")
traj1.add_maneuver(0, maneuver1)
traj1.add_maneuver(1, maneuver2)
t_clock_burn = traj1.get_trajectory_position(0, "last")["t_clock"]
t_clock_arrival = traj1.get_trajectory_position(1, "last")["t_clock"]
traj1.add_trajectory_position(2, t_clock=t_clock_arrival + traj1.orbits[2]["orbit"].T/2, name="End position")

plotter = Plotter(plot3d=True)
plotter.plot_trajectory(traj1, frame="bodycentric", velocities=True, positions=True, orbits=True, maneuvers=True, orbit_arrows=False)

