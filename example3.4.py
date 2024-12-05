from orbitalmechanics import Orbit
import numpy as np
from scipy.optimize import fsolve

"""
A geocentric parabola has a perigee velocity of 10 km/s. How far is the satellite from the center of the 
earth 6 h after perigee passage?
"""

M_earth = 5.974e24 # [kg]
Earth_radius = 6378 # [km]
v = 10 # [km/s]
rp = 2 * Orbit.G * M_earth / v**2 # [km]
t = 6*60*60 # [s]

orbita = Orbit(m1=M_earth, rp=rp, e=1, body1radius=Earth_radius)
orbita.add_position(0, "perigeu")
theta = orbita.theta_from_t(t)
orbita.add_position(theta, "6 horas após perigeu")
r = orbita.r_at_theta(theta)
print(f"Distância do centro da Terra: {r:.4f} km")
orbita.plot(plot_positions=True)


t = orbita.t_from_theta(theta)
print(f"Tempo após perigeu: {t/(60*60):.4f} h")


