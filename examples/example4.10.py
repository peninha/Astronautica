from astronautica import Orbit, Body, Plotter
import numpy as np

"""
Determine the perigee and apogee for an earth satellite whose orbit satisfies all the following conditions:
it is sun synchronous, its argument of perigee is constant, and its period is 3 hours.
"""
earth = Body("earth")
T = 3*60*60 # [s]

#para sun-synchronous, ascending node precisa avançar 360 / 365.2422 deg/day
Omega_dot = 360 / 365.2422 / (24 * 3600) # [deg/s]
Omega_dot_rad = np.radians(Omega_dot)
mu = earth.mu()
a = np.cbrt(mu * T**2 / (4 * np.pi**2))
i = np.degrees(np.arcsin(2/np.sqrt(5)))
i = 180 - i  # Adicionando a solução para o segundo quadrante
e1 = + np.sqrt(1 + np.sqrt(-3/2 * np.sqrt(mu) * earth.J2 * earth.radius**2 * np.cos(np.radians(i)) / (Omega_dot_rad * a**(7/2))))
e2 = - np.sqrt(1 + np.sqrt(-3/2 * np.sqrt(mu) * earth.J2 * earth.radius**2 * np.cos(np.radians(i)) / (Omega_dot_rad * a**(7/2))))
e3 = + np.sqrt(1 - np.sqrt(-3/2 * np.sqrt(mu) * earth.J2 * earth.radius**2 * np.cos(np.radians(i)) / (Omega_dot_rad * a**(7/2))))
e4 = - np.sqrt(1 - np.sqrt(-3/2 * np.sqrt(mu) * earth.J2 * earth.radius**2 * np.cos(np.radians(i)) / (Omega_dot_rad * a**(7/2))))

print(Omega_dot)
print(a)
print(i)
print(f"e1: {e1:.2f}")
print(f"e2: {e2:.2f}")
print(f"e3: {e3:.2f}") #e3 é a solução correta
print(f"e4: {e4:.2f}")

e = e3

orbita = Orbit.from_elements(main_body=earth, a=a, e=e, i=i, Omega0=0, omega0=270, theta0=0)

rp = orbita.rp
ra = orbita.ra

print(f"Perigee: {rp:.2f} km")
print(f"Apogee: {ra:.2f} km")
print(f"Altitude perigee: {earth.altitude(rp):.2f} km")
print(f"Altitude apogee: {earth.altitude(ra):.2f} km")

plotter = Plotter(plot3d=True)
plotter.plot_orbit(orbita, frame="bodycentric")
