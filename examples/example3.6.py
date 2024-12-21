from orbitalmechanics import Orbit, Trajectory
import numpy as np

"""
An earth satellite has an initial true anomaly of θ = 30°, a radius of r0 = 10000 km, and a speed of v0 = 10 km/s.
Use the universal Kepler’s equation to find the change in universal anomaly χ after 1 h and
use that information to determine the true anomaly θ at that time.
"""
M_earth = 5.974e24 # [kg]
Earth_radius = 6378 # [km]

theta0 = 30 # [deg]
t0 = 0 # [s]
t1 = 1*60*60 # [s]
r0 = 10000 # [km]
v0 = 10 # [km/s]

orbita = Orbit.init_from_r_v_theta(m1=M_earth, r=r0, v=v0, theta=theta0, body1radius=Earth_radius)
r_vec_0, v_vec_0 = orbita.state_perifocal(theta0)

trajectory = Trajectory(m1=M_earth, r_vec_0=r_vec_0, v_vec_0=v_vec_0, t0=t0, body1radius=Earth_radius)
Q1 = trajectory.Q_at_delta_t(t1)

print(f"Q1: {Q1}")
#print(f"Posição 60 min após t0: {r_vec_t}")
#print(f"Velocidade 60 min após t0: {v_vec_t}")

#trajectory.plot()
