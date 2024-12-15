from orbitalmechanics import Orbit
import numpy as np

"""
Determine the perigee and apogee for an earth satellite whose orbit satisfies all the following conditions:
it is sun synchronous, its argument of perigee is constant, and its period is 3 hours.
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137
T = 3*60*60 # [s]
J2 = 1.0826e-3
omega_earth = np.degrees(7.292115e-5) # [°/s]

#para sun-synchronous, ascending node precisa avançar 360 / 365.2422 deg/day
Omega_dot = 360 / 365.2422 / (24 * 3600) # [deg/s]
Omega_dot_rad = np.radians(Omega_dot)
mu = Orbit.mu_from_m1_m2(M_earth, 0)
a = np.cbrt(mu * T**2 / (4 * np.pi**2))
i = np.degrees(np.arcsin(2/np.sqrt(5)))
i = 180 - i  # Adicionando a solução para o segundo quadrante
e1 = + np.sqrt(1 + np.sqrt(-3/2 * np.sqrt(mu) * J2 * R_terra**2 * np.cos(np.radians(i)) / (Omega_dot_rad * a**(7/2))))
e2 = - np.sqrt(1 + np.sqrt(-3/2 * np.sqrt(mu) * J2 * R_terra**2 * np.cos(np.radians(i)) / (Omega_dot_rad * a**(7/2))))
e3 = + np.sqrt(1 - np.sqrt(-3/2 * np.sqrt(mu) * J2 * R_terra**2 * np.cos(np.radians(i)) / (Omega_dot_rad * a**(7/2))))
e4 = - np.sqrt(1 - np.sqrt(-3/2 * np.sqrt(mu) * J2 * R_terra**2 * np.cos(np.radians(i)) / (Omega_dot_rad * a**(7/2))))

print(Omega_dot)
print(a)
print(i)
print(f"e1: {e1:.2f}")
print(f"e2: {e2:.2f}")
print(f"e3: {e3:.2f}") #e3 é a solução correta
print(f"e4: {e4:.2f}")

rp = Orbit.rp_from_a_e(a, e3)
ra = Orbit.ra_from_a_e(a, e3)

print(f"Perigee: {rp:.2f} km")
print(f"Apogee: {ra:.2f} km")

orbita = Orbit(ra=ra, rp=rp, i=i, Omega=0, omega=0, theta0=0, m1=M_earth, m2=0)

orbita.add_oblateness_correction(J2, R_terra)
orbita.add_main_body_rotation(omega_earth)

orbita.trajectory(theta1 = 360*20, n_points=21)
orbita.plot(frame="bodycentric",
            plot3d=True, 
            orbit=False,
            points=False,
            positions=True,
            velocities=True,
            trajectory=True)

