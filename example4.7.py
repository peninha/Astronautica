from orbitalmechanics import Orbit
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
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137

h = 80000
e = 1.4
i = 30
Omega = 40
omega = 60
theta = 30

orbita = Orbit(m1=M_earth, m2=0, h=h, e=e, i=i, Omega=Omega, omega=omega, theta=theta, body1radius=R_terra)
r_vec, v_vec = orbita.state_vectors_at_theta(theta, frame="bodycentric")
print(r_vec, v_vec)

orbita.plot(frame="bodycentric", points=True, velocities=True, positions=True, trajectory=False, plot3d=True)

