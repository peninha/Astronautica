from astronautica import Orbit, Body, Plotter
import numpy as np

"""
For a given earth orbit, the elements are:
h = 80000 km2/s
e = 1.4
i = 30°
Ω = 40°
ω = 60°
θ = 30°
Using Algorithm 4.5, find the state vectors r and v in the geocentric equatorial frame.
"""

earth = Body("earth")

h = 80000
e = 1.4
i = 30
Omega = 40
omega = 60
theta0 = 30

orbita = Orbit.from_elements(earth, h=h, e=e, i=i, Omega0=Omega, omega0=omega, theta0=theta0)
r_vec, v_vec = orbita.state_vectors_at_theta(theta0, frame="bodycentric")
print(r_vec, v_vec)

plotter = Plotter(frame="bodycentric", plot3d=True)
plotter.plot_orbit(orbita)
