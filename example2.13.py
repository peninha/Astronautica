from orbitalmechanics import Orbit, Lagrange_mechanics
import numpy as np

"""
An earth satellite moves in the xy plane of an inertial frame with the origin at the earth’s center. Relative to that frame, the
position and velocity of the satellite at time t0 are
r0 = 8182.4i - 6865.9j km
v0 = 0.47572i + 8.8116j km/s
Use Lagrange coefficients to compute the position and velocity vectors after the satellite has traveled through a true anomaly
of 120°
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137 # [km]

# Vetores de posição e velocidade no sistema perifocal
r0_vec = np.array([8182.4, -6865.9, 0])  # km
v0_vec = np.array([0.47572, 8.8116, 0])  # km/s

lagrange = Lagrange_mechanics(m1=M_earth, r_vec_0=r0_vec, v_vec_0=v0_vec, delta_theta=120, body1radius=R_terra)
theta0 = lagrange.theta0
r0 = lagrange.r0
theta = lagrange.theta
r = lagrange.r

lagrange.plot(points=[(r0, theta0), (r, theta)])

print(f"Posição: {lagrange.r_vec}")
print(f"Velocidade: {lagrange.v_vec}")

"""
Find the eccentricity of the orbit in Example 2.13 as well as the true anomaly at the initial time t0 and, hence, the location of
the perigee for this orbit.
"""
print(f"Eccentricidade: {lagrange._e}")
print(f"Anomalia verdadeira no tempo inicial: {theta0}")

orbit = Orbit.init_from_r_vec_v_vec(m1=M_earth, r_vec=r0_vec, v_vec=v0_vec, body1radius=R_terra)
orbit.plot()