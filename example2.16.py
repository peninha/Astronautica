from orbitalmechanics import Three_body_restricted
import numpy as np

"""
Locate the five Lagrange points for the earth–moon system.
m1 = 5.974e24 kg earth
m2 = 7.348e22 kg moon
r12 = 3.844e5 km distance between the earth and moon
"""

M_earth = 5.974e24 # [kg]
M_moon = 7.348e22 # [kg]
r12 = 3.844e5 # [km]

three_body = Three_body_restricted(m1=M_earth, m2=M_moon, r12=r12, body1radius=6378, body2radius=1737)
L1, L2, L3, L4, L5 = three_body.lagrange_points(plot=True)

print(L1)
print(L2)
print(L3)
print(L4)
print(L5)

print(three_body.jacobi_constant(L1, 0))
print(three_body.jacobi_constant(L2, 0))
print(three_body.jacobi_constant(L3, 0))
print(three_body.jacobi_constant(L4, 0))
print(three_body.jacobi_constant(L5, 0))
