from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
A spacecraft returning from a lunar mission approaches earth on a hyperbolic trajectory. At its closest approach A it is at an
altitude of 5000 km, traveling at 10 km/s. At A retrorockets are fired to lower the spacecraft into a 500-km-altitude circular
orbit, where it is to rendezvous with a space station. Find the location of the space station at retrofire so that rendezvous will
occur at B (Fig. 6.6).
"""
earth = Body("spherical_earth")

t0_clock = 0
r1_vec = np.array([-earth.radius_from_altitude(5000), 0, 0])
v1_vec = np.array([0, -10, 0])

orbit_from = Orbit.from_state_vectors(earth, r1_vec, v1_vec, t0_clock=t0_clock)

rp = earth.radius_from_altitude(500)
e = 0.00

orbit_to = Orbit.from_elements(earth, e=e, rp=rp, theta0=0, t0_clock=t0_clock)

maneuver1, maneuver2, delta_v = Maneuver.hohmann_transfer(orbit_from, orbit_to, theta_burn=0)

delta_v1 = maneuver1.delta_v
delta_v2 = maneuver2.delta_v
print(f"Delta-v 1: {delta_v1} km/s")    
print(f"Delta-v 2: {delta_v2} km/s")
print(f"Total delta-v: {delta_v} km/s")

traj1 = Trajectory(orbit_from, orbit_name="Spacecraft orbit", position_name="Spacecraft burn position")
traj1.add_trajectory_position(0, t_clock=t0_clock-3000, name="Spacecraft initial position")
traj1.add_maneuver(0, maneuver1, position_name="Spacecraft burn position")
traj1.add_maneuver(1, maneuver2, position_name="Spacecraft rendezvous position")


t_clock_arrival = traj1.get_trajectory_position(1, position_index="last")['t_clock']
station_orbit = Orbit.from_elements(earth, e=0, rp=rp, theta0=0, t0_clock=t_clock_arrival)
print(f"t_clock_arrival: {t_clock_arrival}")
print(f"station orbit period: {station_orbit.T}")

traj2 = Trajectory(station_orbit, orbit_name="Station orbit", position_name="Station rendezvous position")
traj2.add_trajectory_position(0, t_clock=t0_clock, name="Station initial position")
station_theta0 = station_orbit.theta_at_t_clock(t0_clock)
print(f"station_theta0: {station_theta0}")


plotter = Plotter(plot3d=False)
plotter.plot_trajectories([traj1, traj2],
                        frame="bodycentric",
                        points=True,
                        velocities=True,
                        positions=True,
                        orbits=False,
                        maneuvers=True,
                        v_scale=1)

