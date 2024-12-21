from orbitalmechanics import Orbit
import numpy as np

"""
The position of an earth satellite is first determined to be
r1 = (5000, 10000, 2100) km
After 1 h the position vector is
r2 = (-14600, 2500, 7000) km
Determine the orbital elements and find the perigee altitude and the time since perigee passage of the first sighting.
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137

r1 = np.array([5000, 10000, 2100])
r2 = np.array([-14600, 2500, 7000])
delta_t = 1 * 3600 # [s]
orbita = Orbit.init_from_2_vectors_and_delta_time(r1, r2, delta_t, m1=M_earth, m2=0, body1radius=R_terra)

print(orbita)

# calculate perigee
rp = orbita.rp 

# calculate time since perigee passage
first_position = min(position['theta'] for position in orbita.positions)
t_orbit_first_position = orbita.t_orbit_at_theta(first_position)

print(f"Perigee altitude: {rp - R_terra} km")
print(f"Time since perigee passage: {t_orbit_first_position} s")

orbita.plot(frame="bodycentric", points=True, velocities=True, positions=True, trajectory=False, plot3d=True)