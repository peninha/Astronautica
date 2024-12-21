from orbitalmechanics import Orbit
import numpy as np

"""
A geocentric elliptical orbit has a perigee radius of 9600 km and an apogee radius of 21,000 km, as shown in Fig. 3.8.
Calculate the time to fly from perigee P to a true anomaly of 120°.
"""

M_earth = 5.974e24 # [kg]
Earth_radius = 6378 # [km]
rp = 9600 # [km]
ra = 21000 # [km]
theta0 = 0 # [°]
theta1 = 120 # [°]

orbita = Orbit(m1=M_earth, rp=rp, ra=ra, body1radius=Earth_radius)
orbita.add_orbital_position(name="sonda0", theta=theta0)
orbita.add_orbital_position(name="sonda1", theta=theta1)
t = orbita.t_at_theta(theta1) - orbita.t_at_theta(theta0)
print(f"Time to fly from perigee to a true anomaly of {theta1}°: {t/3600:.2f} hours")
orbita.plot(plot_positions=True)