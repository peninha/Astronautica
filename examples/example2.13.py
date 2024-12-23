from astronautica import Orbit, Body, Plotter
import numpy as np

"""
An earth satellite moves in the xy plane of an inertial frame with the origin at the earth’s center. Relative to that frame, the
position and velocity of the satellite at time t0 are
r0 = 8182.4i - 6865.9j km
v0 = 0.47572i + 8.8116j km/s
Use Lagrange coefficients to compute the position and velocity vectors after the satellite has traveled through a true anomaly
of 120°
"""
earth = Body(name="earth")

# Vetores de posição e velocidade no sistema bodycentric
r0_vec_bc = np.array([8182.4, -6865.9, 0])  # km
v0_vec_bc = np.array([0.47572, 8.8116, 0])  # km/s
delta_theta = 120 # [°]
r0 = np.linalg.norm(r0_vec_bc)

orbita = Orbit.from_state_vectors(earth, r_vec_bc=r0_vec_bc, v_vec_bc=v0_vec_bc)

theta0 = orbita.orbital_positions[0]['theta']
orbita.add_orbital_position(theta=theta0+delta_theta, name="Position 1")
r1_vec_bc, v1_vec_bc = orbita.state_vectors_at_theta(theta0+delta_theta, "bodycentric")

plot = Plotter(frame="bodycentric", plot3d=False)
plot.plot_orbit(orbit=orbita)

print(f"Posição0: {r0_vec_bc}")
print(f"Velocidade0: {v0_vec_bc}")
print(f"Posição1: {r1_vec_bc}")
print(f"Velocidade1: {v1_vec_bc}")

"""
Find the eccentricity of the orbit in Example 2.13 as well as the true anomaly at the initial time t0 and, hence, the location of
the perigee for this orbit.
"""
print(f"Eccentricidade: {orbita.e}")
print(f"Anomalia verdadeira no tempo inicial: {theta0}")