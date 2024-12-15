from orbitalmechanics import Orbit
import numpy as np

"""
Molniya orbit
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137
T = 12*60*60 - 130 # [s]
J2 = 1.0826e-3
omega_earth = np.degrees(7.292115e-5) # [deg/s]

mu = Orbit.mu_from_m1_m2(M_earth, 0)
i = np.degrees(np.arcsin(2/np.sqrt(5)))
e = 0.72
a = np.cbrt(mu * T**2 / (4 * np.pi**2))
Omega_dot = np.degrees(-3/2 * np.sqrt(mu) * J2 * R_terra**2 * np.cos(np.radians(i)) / (a**(7/2) * (1 - e**2)**2))

print('Omega_dot: ', Omega_dot)
print('omega_earth: ', omega_earth)
print('a: ', a)
print('i: ', i)
print('e: ', e)

rp = Orbit.rp_from_a_e(a, e)
ra = Orbit.ra_from_a_e(a, e)

print(f"Perigee: {rp - R_terra:.2f} km")
print(f"Apogee: {ra - R_terra:.2f} km")

orbita = Orbit(ra=ra, rp=rp, i=i, Omega=0, omega=270, theta0=0, m1=M_earth, m2=0)
orbita.add_oblateness_correction(J2, R_terra)
orbita.add_main_body_rotation(omega_earth)

orbita.trajectory(theta1 = 360*50, n_points=51)
orbita.plot(frame="rotatingBodycentric",
            plot3d=True, 
            orbit=False,
            points=False,
            positions=True,
            velocities=True,
            trajectory=True)

