from astronautica import Orbit, Body, Plotter
import numpy as np

"""
Given the state vector,
r = - 6045I - 3490J + 2500K km
v = - 3.457I + 6.618J + 2.533K km/s
find the orbital elements h, i, Ω, e, ω, and θ.
"""

earth = Body("earth")

r_vec = np.array([-6045, -3490, 2500]) # [km]
v_vec = np.array([-3.457, 6.618, 2.533]) # [km/s]

orbita = Orbit.from_state_vectors(earth, r_vec, v_vec)

print(f"h: {orbita.h:.4f} km²/s")
print(f"i: {orbita.i:.4f}°")
print(f"Ω: {orbita.Omega0:.4f}°")
print(f"e: {orbita.e:.4f}")
print(f"ω: {orbita.omega0:.4f}°")
print(f"θ: {orbita.theta0:.4f}°")


theta1 = 180
orbita.add_orbital_position(name="Apoapsis", theta=theta1)

plotter = Plotter(frame="bodycentric", plot3d=True)
plotter.plot_orbit(orbita)
