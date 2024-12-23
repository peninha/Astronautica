from astronautica import Orbit, Body, Plotter
import numpy as np
"""
At time t0, the state vector of an earth satellite is
r0 = 1600I + 5310J + 3800K km
v0 = -7.350I + 0.4600J + 2.470K km/s
Determine the position and velocity vectors 3200s later and plot
the orbit in three dimensions.
"""

earth = Body("earth")

r_vec_0 = np.array([1600, 5310, 3800]) # [km]
v_vec_0 = np.array([-7.350, 0.4600, 2.470]) # [km/s]

t0 = 0 # [s]
t1 = 3200 # [s]
orbita = Orbit.from_state_vectors(earth, r_vec_0, v_vec_0, t0_clock=t0)
orbita.add_orbital_position(name="Position 1", t_clock=t1)
r_vec_1, v_vec_1 = orbita.state_vectors_at_t_clock(t1, "bodycentric")

print(f"Posição 3200s após t0: {r_vec_1}")
print(f"Velocidade 3200s após t0: {v_vec_1}")

plotter = Plotter(plot3d=True)
plotter.plot_orbit(orbita, frame="bodycentric")
