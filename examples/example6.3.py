from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
Find the total delta-v requirement for bielliptic Hohmann transfer from a geocentric circular orbit of 7000 km
radius to one of 105,000 km radius. Let the apogee of the first ellipse be 210,000 km. Compare the delta-v
schedule and total flight time with that for an ordinary single Hohmann transfer ellipse (see Fig. 6.9).
"""
earth = Body("spherical_earth")

t0_clock = 0
rp1 = 7000
e1 = 0
orbit_from = Orbit.from_elements(earth, e=e1, rp=rp1, theta0=0, t0_clock=t0_clock)

rp2 = 105000
e2 = 0

orbit_to = Orbit.from_elements(earth, e=e2, rp=rp2, theta0=0, t0_clock=t0_clock)

maneuver1, maneuver2, maneuver3, delta_v = Maneuver.bielliptic_transfer(orbit_from, orbit_to, apoapsis=210000, theta_burn=0)

delta_v1 = maneuver1.delta_v
delta_v2 = maneuver2.delta_v
delta_v3 = maneuver3.delta_v
print(f"Delta-v 1: {delta_v1} km/s")    
print(f"Delta-v 2: {delta_v2} km/s")
print(f"Delta-v 3: {delta_v3} km/s")
print(f"Total delta-v: {delta_v} km/s")

traj1 = Trajectory(orbit_from, orbit_name="Spacecraft orbit", position_name="Spacecraft burn position")
traj1.add_trajectory_position(0, t_clock=t0_clock-3000, name="Spacecraft initial position")
traj1.add_maneuver(0, maneuver1, position_name="Spacecraft burn position")
traj1.add_maneuver(1, maneuver2, position_name="Spacecraft transfer apoapsis position")
traj1.add_maneuver(2, maneuver3, position_name="Spacecraft arrival position")
traj1.add_trajectory_position(3, t_clock=traj1.get_trajectory_position(2, position_index="last")['t_clock']+10000, name="Spacecraft final position")

plotter = Plotter(plot3d=False)
plotter.plot_trajectory(traj1,
                        frame="bodycentric",
                        points=True,
                        velocities=True,
                        positions=True,
                        orbits=False,
                        maneuvers=True,
                        time_step=700,
                        v_scale=1)

