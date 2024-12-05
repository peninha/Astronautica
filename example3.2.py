from orbitalmechanics import Orbit
import numpy as np

"""
A geocentric elliptical orbit has a perigee radius of 9600 km and an apogee radius of 21,000 km, as shown in Fig. 3.8.
Find the true anomaly at 3 hours after perigee passage.
"""

M_earth = 5.974e24 # [kg]
Earth_radius = 6378 # [km]
rp = 9600 # [km]
ra = 21000 # [km]
theta0 = 0 # [°]
t = 3*3600 # [s]

orbita = Orbit(m1=M_earth, rp=rp, ra=ra, body1radius=Earth_radius)
orbita.add_position(name="sonda0", theta=theta0)
theta1 = orbita.theta_from_t(t)
orbita.add_position(name="sonda1", theta=theta1)
print(f"True anomaly at {t/3600:.2f} hours after perigee passage: {theta1:.2f}°")
orbita.plot(plot_positions=True)