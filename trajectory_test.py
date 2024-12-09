from orbitalmechanics import Orbit
import numpy as np

M_earth = 5.974e24 # [kg]
Earth_radius = 6378 # [km]

e = 1.2
rp = 7000 # [km]

t0_clock = 120 # [s]
t1_clock = 36000 # [s]
theta0 = -130 # [rad]
theta1 = 129 # [rad]
n_points = 60

orbita = Orbit(rp=rp, e=e, m1=M_earth, body1radius=Earth_radius)
orbita.add_zero_state_from_theta(theta0=theta0, t0_clock=t0_clock)
orbita.trajectory(theta1=theta1, n_points=n_points)
orbita.plot(orbit=False)

print(orbita)
