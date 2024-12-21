from orbitalmechanics import Three_body_restricted
import numpy as np

"""
The earth-orbiting spacecraft has a relative burnout velocity vbo at an altitude of d = 200 km
on a radial for which ϕ = -90° relative to the Earth. Find the value of vbo for each Lagrange point.
"""

M_earth = 5.974e24 # [kg]
M_moon = 7.348e22 # [kg]
r12 = 3.844e5 # [km]
d = 200 # [km]
Earth_radius = 6378 # [km]
Moon_radius = 1737 # [km]

three_body = Three_body_restricted(m1=M_earth, m2=M_moon, r12=r12, body1radius=Earth_radius, body2radius=Moon_radius)
L1, L2, L3, L4, L5 = three_body.lagrange_points(add_points=True)
three_body.plot(frame="rotatingBarycentric")

C1 = three_body.jacobi_constant(L1, 0)
C2 = three_body.jacobi_constant(L2, 0)
C3 = three_body.jacobi_constant(L3, 0)
C4 = three_body.jacobi_constant(L4, 0)
C5 = three_body.jacobi_constant(L5, 0)

r_vec = np.array([-three_body.pi2*three_body.r12, -(d + Earth_radius), 0])
print(r_vec)

v1 = three_body.v_for_C(r_vec, C1)
v2 = three_body.v_for_C(r_vec, C2)
v3 = three_body.v_for_C(r_vec, C3)
v4 = three_body.v_for_C(r_vec, C4)
v5 = three_body.v_for_C(r_vec, C5)

print(f"C1 ({C1:.5f}): vbo = {v1:.5f} km/s")
print(f"C2 ({C2:.5f}): vbo = {v2:.5f} km/s")
print(f"C3 ({C3:.5f}): vbo = {v3:.5f} km/s")
print(f"C4 ({C4:.5f}): vbo = {v4:.5f} km/s")
print(f"C5 ({C5:.5f}): vbo = {v5:.5f} km/s")
