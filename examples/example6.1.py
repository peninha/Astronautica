from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
A 2000-kg spacecraft is in a 480 km by 800 km earth orbit (orbit 1 in Fig. 6.3). Find
(a) The Δv required at perigee A to place the spacecraft in a 480 km by 16,000 km transfer ellipse (orbit 2).
(b) The Δv (apogee kick) required at B of the transfer orbit to establish a circular orbit of 16,000 km altitude (orbit 3).
(c) The total required propellant if the specific impulse is 300 s
"""
earth = Body("spherical_earth")

rp = earth.radius_from_altitude(480)
ra = earth.radius_from_altitude(800)
theta0 = -180
i = 0
Omega0 = 0
omega0 = 15

orbit_from = Orbit.from_apsis(earth, rp, ra, i=i, theta0=theta0, Omega0=Omega0, omega0=omega0)

rp = earth.radius_from_altitude(16000)
e = 0.00

orbit_to = Orbit.from_elements(earth, e=e, rp=rp, i=i, theta0=0, Omega0=Omega0, omega0=omega0)

maneuver1, maneuver2, delta_v = Maneuver.hohmann_transfer(orbit_from, orbit_to, theta_burn=0)

delta_v1 = maneuver1.delta_v
delta_v2 = maneuver2.delta_v
print(f"Delta-v 1: {delta_v1} km/s")    
print(f"Delta-v 2: {delta_v2} km/s")
print(f"Total delta-v: {delta_v} km/s")

traj1 = Trajectory(orbit_from)
traj1.add_maneuver(0, maneuver1)
traj1.add_maneuver(1, maneuver2)
traj1.add_trajectory_position(2, t_clock=traj1.get_trajectory_position(2, position_index="last")['t_clock'] + 6000, name="Final position")

traj2 = Trajectory(orbit_to)
traj2.add_trajectory_position(0, t_clock=3000, name="test")

plotter = Plotter(plot3d=False)
plotter.plot_trajectories([traj1, traj2],
                        frame="bodycentric",
                        points=True,
                        velocities=True,
                        positions=True,
                        orbits=True,
                        maneuvers=True,
                        v_scale=5)

