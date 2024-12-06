from orbitalmechanics import Trajectory
import numpy as np

"""
At time t0, the state vector of an earth satellite is
r0 = 1600I + 5310J + 3800K km
v0 = -7.350I + 0.4600J + 2.470K km/s
Determine the position and velocity vectors 3200 s later and plot
the orbit in three dimensions.
"""

M_earth = 5.974e24 # [kg]
Earth_radius = 6378 # [km]

r_vec_0 = np.array([1600, 5310, 3800]) # [km]
v_vec_0 = np.array([-7.350, 0.4600, 2.470]) # [km/s]


t = 3200 # [s]
trajectory = Trajectory(m1=M_earth, r_vec_0=r_vec_0, v_vec_0=v_vec_0, t0=t, body1radius=Earth_radius)
r_vec_t, v_vec_t = trajectory.state_at_Q(trajectory.Q_at_delta_t(t), t)

print(f"Posição 3200 s após t0: {r_vec_t}")
print(f"Velocidade 3200 s após t0: {v_vec_t}")

#trajectory.plot()
