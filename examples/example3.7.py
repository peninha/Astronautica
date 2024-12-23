from astronautica import Orbit, Body, Plotter
import numpy as np

import orbitalmechanics

"""
An earth satellite moves in the xy plane of an inertial frame with origin at the earth’s center.
Relative to that frame, the position and velocity of the satellite at time t0 are
r0 =  7000i - 12124j (km)
v0 = 2.6679i + 4.6210j (km/s)
Compute the position and velocity vectors of the satellite 60 min later.
"""

earth = Body(name="earth")
t0 = 0 # [s]
t1 = 60*60 # [s]

r_vec_0 = np.array([7000, -12124, 0]) # [km]
v_vec_0 = np.array([2.6679, 4.6210, 0]) # [km/s]

orbita = Orbit.from_state_vectors(earth, r_vec_0, v_vec_0, t0_clock=t0)
orbita.add_orbital_position(t_clock=t1, name="60 min após t0")
r_vec_1, v_vec_1 = orbita.state_vectors_at_t_clock(t1, frame="bodycentric")

print(f"Posição 60 min após t0: {r_vec_1}")
print(f"Velocidade 60 min após t0: {v_vec_1}")

plot = Plotter(frame="bodycentric", plot3d=False)
plot.plot_orbit(orbit=orbita)
