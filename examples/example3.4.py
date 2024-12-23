from astronautica import Orbit, Body, Plotter
import numpy as np
from scipy.optimize import fsolve

"""
A geocentric parabola has a perigee velocity of 10 km/s. How far is the satellite from the center of the 
earth 6 h after perigee passage?
"""

earth = Body(name="earth")
e = 1
v = 10 # [km/s]
rp = 2 * Orbit.G * earth.mass / v**2 # [km]
t = 6*60*60 # [s]

orbita = Orbit.from_elements(earth, rp=rp, e=e, theta0=0)
theta = orbita.theta_at_t_clock(t)
orbita.add_orbital_position(t, name="6 horas após perigeu")
r = orbita.r_at_theta(theta)
print(f"Distância do centro da Terra: {r:.4f} km")

plot = Plotter(frame="bodycentric", plot3d=True)
plot.plot_orbit(orbit=orbita)


t = orbita.t_clock_at_theta(theta)
print(f"Tempo após perigeu: {t/(60*60):.4f} h")


