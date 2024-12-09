from orbitalmechanics import Orbit
import numpy as np
"""
At time t0, the state vector of an earth satellite is
r0 = 1600I + 5310J + 3800K km
v0 = -7.350I + 0.4600J + 2.470K km/s
Determine the position and velocity vectors 3200s later and plot
the orbit in three dimensions.
"""

M_earth = 5.974e24 # [kg]
Earth_radius = 6378 # [km]

r_vec_0 = np.array([1600, 5310, 3800]) # [km]
v_vec_0 = np.array([-7.350, 0.4600, 2.470]) # [km/s]

t0 = 0 # [s]
t1 = 3200 # [s]
orbita = Orbit.init_from_state_vectors(r_vec_0, v_vec_0, m1=M_earth, body1radius=Earth_radius, starting_point=True, t0=t0)
r_vec_1, v_vec_1 = orbita.state_at_t(t1)

orbita.trajectory(t1_clock=t1)
orbita.plot(orbit=False)

print(f"Posição 3200s após t0: {r_vec_1}")
print(f"Velocidade 3200s após t0: {v_vec_1}")
