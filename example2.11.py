from orbitalmechanics import Orbit
import numpy as np

"""
An earth orbit has an eccentricity of 0.3, an angular momentum of 60000 km2/s, 
and a true anomaly of 120°. What are the position vector r and velocity vector v
in the perifocal frame of reference?
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137 # [km]

orbit = Orbit(m1=M_earth, e=0.3, h=60000, body1radius=R_terra)

print(orbit.r_vec_at_theta(120))
print(orbit.v_vec_at_theta(120))

r = orbit.r_at_theta(120)
orbit.add_orbital_position(120, name='r')
orbit.plot(plot_velocities=True)