from astronautica import Orbit, Body, Plotter
import numpy as np

"""
A meteoroid is sighted at an altitude of 267,000 km. After 13.5 h and a change in true anomaly of 5°,
the altitude is observed to be 140,000 km.
Calculate the perigee altitude and the time to perigee after the second sighting.
"""
earth = Body("earth")

r1 = earth.radius_from_altitude(267000)
r2 = earth.radius_from_altitude(140000)
delta_t = 13.5 * 3600 # [s]
delta_theta = 5 # [deg]

orbita = Orbit.from_2_radii_delta_t_delta_theta(earth, r1, r2, delta_t, delta_theta)

#print(orbita)

rp = orbita.rp

print(f"Perigee altitude: {earth.altitude(rp)} km")

# calculate time to perigee after the second sighting
second_position = orbita.orbital_positions[0]['theta']
t_orbit_second_position = orbita.t_orbit_at_theta(second_position)

print(f"Time to perigee after the second sighting: {t_orbit_second_position} s")

plot = Plotter(plot3d=True)
plot.plot_orbit(orbita, frame="bodycentric", points=True, velocities=True, positions=True)
