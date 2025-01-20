from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
Find the delta-v required to transfer a satellite from a circular,
300-km-altitude low earth orbit of 28° inclination to a geostationary equatorial orbit.
Circularize and change the inclination at altitude. Compare that delta-v requirement with
the one in which the plane change is done in the low earth orbit.
"""
earth = Body("spherical_earth")

e = 0
rp1 = earth.radius_from_altitude(300)
i1 = 28
rp2 = earth.geosynchronous_radius()
i2 = 0

orbit_from = Orbit.from_elements(earth, e=e, rp=rp1, i=i1, theta0=0, t0_clock=0, name="Orbit from")
orbit_to = Orbit.from_elements(earth, e=e, rp=rp2, i=i2, theta0=0, t0_clock=0, name="Orbit to")

theta_burn = 0.001
theta_arrival = 180

maneuver1, maneuver2, delta_t, delta_v = Maneuver.orbit_pos_to_orbit_pos_maneuver(orbit_from, orbit_to, theta_burn, theta_arrival)

traj1 = Trajectory(orbit_from)
traj1.add_trajectory_position(0, -orbit_from.T/2)
traj1.add_maneuver(0, maneuver1)
traj1.add_maneuver(1, maneuver2)
traj1.add_trajectory_position(2, delta_t + 5000)

print(f"Delta-v required: {delta_v:.5f} km/s")

plotter = Plotter(plot3d=True)
plotter.plot_trajectory(traj1, orbits=True)