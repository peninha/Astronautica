from orbitalmechanics import Three_body_restricted
import numpy as np

"""
Plot lagrandian points of the Sun-Earth system.
"""

M_earth = 5.974e24 # [kg]
M_sun = 1.989e30 # [kg]
r12 = 1.496e8 # [km]
Earth_radius = 6378 # [km]
Sun_radius = 696340 # [km]

three_body = Three_body_restricted(m1=M_sun, m2=M_earth, r12=r12, body1radius=Sun_radius, body2radius=Earth_radius)
three_body.lagrange_points(add_points=True)

three_body.plot(frame="rotatingBarycentric")