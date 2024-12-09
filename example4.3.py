from orbitalmechanics import Orbit
import numpy as np

"""
Given the state vector,
r = - 6045I - 3490J + 2500K km
v = - 3.457I + 6.618J + 2.533K km/s
find the orbital elements h, i, Ω, e, ω, and θ.
"""

M_earth = 5.974e24 # [kg]
Earth_radius = 6378 # [km]

r_vec = np.array([-6045, -3490, 2500]) # [km]
v_vec = np.array([-3.457, 6.618, 2.533]) # [km/s]

orbita = Orbit.init_from_state_vectors(r_vec, v_vec, m1=M_earth, body1radius=Earth_radius)

print(f"h: {orbita.h:.4f} km²/s")
print(f"i: {orbita.i:.4f}°")
print(f"Ω: {orbita.Omega:.4f}°")
print(f"e: {orbita.e:.4f}")
print(f"ω: {orbita.omega:.4f}°")
print(f"θ: {orbita.theta:.4f}°")

orbita.plot()
