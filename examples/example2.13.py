from orbitalmechanics import Orbit
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
delta_theta = 120 # [°]
r0 = np.linalg.norm(r0_vec)

orbita = Orbit.init_from_state_vectors(m1=M_earth, r_vec=r0_vec, v_vec=v0_vec, body1radius=R_terra)
theta0 = orbita.theta_at_state_vectors(r0_vec, v0_vec)
r_vec, v_vec = orbita.state_after_delta_theta(r0_vec, v0_vec, delta_theta)
theta1 = orbita.theta_at_state_vectors(r_vec, v_vec)
orbita.add_orbital_position(theta1, name="theta = 120°")
orbita.plot()

print(f"Posição: {r_vec}")
print(f"Velocidade: {v_vec}")

"""
Find the eccentricity of the orbit in Example 2.13 as well as the true anomaly at the initial time t0 and, hence, the location of
the perigee for this orbit.
"""
print(f"Eccentricidade: {orbita.e}")
print(f"Anomalia verdadeira no tempo inicial: {theta0}")