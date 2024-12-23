from astronautica import Orbit, Body, Plotter
import numpy as np

"""
A geocentric elliptical orbit has a perigee radius of 9600 km and an apogee radius of 21,000 km, as shown in Fig. 3.8.
Calculate the time to fly from perigee P to a true anomaly of 120°.
"""

earth = Body(name="earth")
rp = 9600 # [km]
ra = 21000 # [km]
theta0 = 0 # [°]
theta1 = 120 # [°]

orbita = Orbit.from_elements(earth, rp=rp, ra=ra, theta0=theta0)
orbita.add_orbital_position(name="Position 1", theta=theta1)
t = orbita.orbital_positions[1]['t_clock'] - orbita.orbital_positions[0]['t_clock']
print(f"Time to fly from perigee to a true anomaly of {theta1}°: {t/3600:.2f} hours")
plot = Plotter(frame="bodycentric", plot3d=True)
plot.plot_orbit(orbit=orbita)