from astronautica import Orbit, Body, Plotter
import numpy as np

"""
A geocentric elliptical orbit has a perigee radius of 9600 km and an apogee radius of 21,000 km, as shown in Fig. 3.8.
Find the true anomaly at 3 hours after perigee passage.
"""

earth = Body(name="earth")
rp = 9600 # [km]
ra = 21000 # [km]
theta0 = 0 # [°]
t = 3*3600 # [s]

orbita = Orbit.from_elements(earth, rp=rp, ra=ra, theta0=theta0)
orbita.add_orbital_position(name="Position 1", t_clock=t)
theta1 = orbita.orbital_positions[1]['theta']
print(f"True anomaly at {t/3600:.2f} hours after perigee passage: {theta1:.2f}°")

plot = Plotter(frame="bodycentric", plot3d=True)
plot.plot_orbit(orbit=orbita)