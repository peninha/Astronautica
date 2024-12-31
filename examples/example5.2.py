from astronautica import Orbit, Body, Plotter
import numpy as np

"""
The position of an earth satellite is first determined to be
r1 = (5000, 10000, 2100) km
After 1 h the position vector is
r2 = (-14600, 2500, 7000) km
Determine the orbital elements and find the perigee altitude and the time since perigee passage of the first sighting.
"""

earth = Body("earth")

r1 = np.array([5000, 10000, 2100])
r2 = np.array([-14600, 2500, 7000])
delta_t = 1 * 3600 # [s]
orbita = Orbit.from_2_vectors_and_delta_time(earth, r1, r2, delta_t)

#print(orbita)

# calculate perigee
rp = orbita.rp 

# calculate time since perigee passage
first_position = min(position['theta'] for position in orbita.orbital_positions)
t_orbit_first_position = orbita.t_orbit_at_theta(first_position)

print(f"Perigee altitude: {earth.altitude(rp)} km")
print(f"Time since perigee passage: {t_orbit_first_position} s")

plot = Plotter(plot3d=True)
plot.plot_orbit(orbita, frame="bodycentric", points=True, velocities=True, positions=True)