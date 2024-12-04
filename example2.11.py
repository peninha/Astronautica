from orbitalmechanics import Orbit
import numpy as np

"""
An earth orbit has an eccentricity of 0.3, an angular momentum of 60, 000 km2/s, 
and a true anomaly of 120°. What are theposition vector r and velocity vector v
in the perifocal frame of reference?
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137 # [km]

orbit = Orbit(m1=M_earth, e=0.3, h=60000, body1radius=R_terra)

print(orbit.r_vec_perifocal(120))
print(orbit.v_vec_perifocal(120))

r = orbit.r_at_theta(120)

orbit.plot(points=[(r, 120)])