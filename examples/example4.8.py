from astronautica import Orbit, Body, Plotter, Trajectory
import numpy as np

"""
A spacecraft is in a 280 km by 400 km orbit with an inclination of 51.43°.
Find the rates of node regression and perigee advance.
"""

earth = Body("earth")
earth.J2 = earth.J2 * 1
rp = earth.radius_from_altitude(280)
ra = earth.radius_from_altitude(400)
i = 51.43

theta0 = 30
theta1 = 50*360 + theta0
t0_clock = 400

orbita = Orbit.from_elements(earth, rp=rp, ra=ra, i=i, Omega0=40, omega0=90, theta0=theta0)

print(orbita.Omega_dot*24*60*60)
print(orbita.omega_dot*24*60*60)

plotter = Plotter(plot3d=True)
#plotter.plot_orbit(orbita, frame="bodycentric")

trajetoria = Trajectory(orbita)
trajetoria.add_trajectory_position(0, t_clock=t0_clock+25000, name="Initial Position")
plotter.plot_trajectory(trajetoria, time_step=60, frame="bodycentric", orbits=True, points=True, velocities=True, positions=True, maneuvers=False, groundtrack=False)

