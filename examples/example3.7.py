from orbitalmechanics import Trajectory
import numpy as np

"""
An earth satellite moves in the xy plane of an inertial frame with origin at the earth’s center.
Relative to that frame, the position and velocity of the satellite at time t0 are
r0 =  7000i - 12124j (km)
v0 = 26679i + 46210j (km/s)
Compute the position and velocity vectors of the satellite 60 min later.
"""

M_earth = 5.974e24 # [kg]
Earth_radius = 6378 # [km]
t0 = 0 # [s]
t1 = 60*60 # [s]

r_vec_0 = np.array([7000, -12124]) # [km]
v_vec_0 = np.array([26679, 46210]) # [km/s]

trajectory = Trajectory(m1=M_earth, r_vec_0=r_vec_0, v_vec_0=v_vec_0, t0=t0, body1radius=Earth_radius)
r_vec_t, v_vec_t = trajectory.state_at_time(t1)

print(f"Posição 60 min após t0: {r_vec_t}")
print(f"Velocidade 60 min após t0: {v_vec_t}")

#trajectory.plot()
