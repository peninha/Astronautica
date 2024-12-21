from orbitalmechanics import Orbit
import numpy as np

"""
A satellite is to be launched into a sun-synchronous circular orbit with a 
period of 100 min. Determine the required altitude and inclination of its orbit.
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137

J2 = 1.0826e-3

T = 100 * 60 # [s]

delta_theta_day = 360 / 365.2422 # [deg/day]
Omega_dot = delta_theta_day / (24 * 3600) # [deg/s]
Omega_dot_rad = np.radians(Omega_dot)
mu = Orbit.mu_from_m1_m2(M_earth, 0)
r = np.cbrt(mu * T**2 / (4 * np.pi**2))
e = 0
a = r

print(f"Altitude: {r - R_terra:.2f} km")

i = np.degrees(np.arccos(-2/3 * Omega_dot_rad * a**(7/2) * (1 - e**2)**2 / (np.sqrt(mu) * J2 * R_terra**2)))

print(f"Inclination: {i:.2f} deg")

orbita = Orbit(a=a, e=e, i=i, Omega=0, omega=0, theta0=0, m1=M_earth, m2=0)
orbita.add_oblateness_correction(J2, R_terra)

orbita.trajectory(theta1 = 360*5, n_points=500)
orbita.plot(frame="bodycentric",
            plot3d=True, 
            orbit=False,
            points=False,
            positions=True,
            velocities=True,
            trajectory=True)




