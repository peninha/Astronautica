from astronautica import Orbit, Body, Plotter
import numpy as np

"""
A satellite is to be launched into a sun-synchronous circular orbit with a 
period of 100 min. Determine the required altitude and inclination of its orbit.
"""
earth = Body("earth")
T = 100 * 60 # [s]

delta_theta_day = 360 / 365.2422 # [deg/day]
Omega_dot = delta_theta_day / (24 * 3600) # [deg/s]
Omega_dot_rad = np.radians(Omega_dot)
mu = earth.mu()
r = np.cbrt(mu * T**2 / (4 * np.pi**2))
e = 0
a = r

print(f"Altitude: {earth.altitude(r):.2f} km")

i = np.degrees(np.arccos(-2/3 * Omega_dot_rad * a**(7/2) * (1 - e**2)**2 / (np.sqrt(mu) * earth.J2 * earth.radius**2)))

print(f"Inclination: {i:.2f} deg")

orbita = Orbit.from_elements(main_body=earth, a=a, e=e, i=i, Omega0=0, omega0=0, theta0=0)

plotter = Plotter(plot3d=True)
plotter.plot_orbit(orbita)




