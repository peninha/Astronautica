from orbitalmechanics import Orbit
import numpy as np

"""
A spacecraft is in a 280 km by 400 km orbit with an inclination of 51.43°. Find the rates of node regression and perigee
advance.
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137

J2 = 1.0826e-3
#J2 = 10*J2

rp = 280 + R_terra
ra = 400 + R_terra
i = 51.43

theta0 = 30
theta1 = 50*360 + theta0
t0_clock = 400

orbita = Orbit(m1=M_earth, m2=0, rp=rp, ra=ra, Omega=40, i=i, omega=90, theta0=theta0, t0_clock=t0_clock, body1radius=R_terra)
orbita.add_oblateness_correction(J2, body1radius=R_terra)

orbita.trajectory(theta1 = theta1, n_points=101)

orbita.plot(frame="bodycentric",
            plot3d=True, 
            orbit=False,
            points=True,
            positions=True,
            velocities=True,
            trajectory=True)

print(orbita.Omega_dot*24*60*60)
print(orbita.omega_dot*24*60*60)