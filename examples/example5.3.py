from orbitalmechanics import Orbit
import numpy as np

"""
A meteoroid is sighted at an altitude of 267,000 km. After 13.5 h and a change in true anomaly of 5°,
the altitude is observed to be 140,000 km.
Calculate the perigee altitude and the time to perigee after the second sighting.
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137

r1 = 267000 + R_terra
r2 = 140000 + R_terra
delta_t = 13.5 * 3600 # [s]
delta_theta = 5 # [deg]

orbita = Orbit.init_from_2_radii_delta_t_delta_theta(r1, r2, delta_t, delta_theta, m1=M_earth, m2=0, body1radius=R_terra)

print(orbita)

rp = orbita.rp

print(f"Perigee altitude: {rp - R_terra} km")

# calculate time to perigee after the second sighting
second_position = orbita.positions[0]['theta']
t_orbit_second_position = orbita.t_orbit_at_theta(second_position)

print(f"Time to perigee after the second sighting: {t_orbit_second_position} s")

orbita.trajectory(theta1=-10, n_points=20)
orbita.plot(frame='perifocal_t0', plot3d=True, groundtrack=False, orbit=True)

