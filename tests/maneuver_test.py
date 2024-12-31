from astronautica import Orbit, Body, Plotter, Frame, Maneuver, Trajectory
import numpy as np

earth = Body(name="earth")

rp = 7000
ra = 35000
theta0 = -150
t0_clock = 0

orbita = Orbit.from_elements(earth, rp=rp, ra=ra, Omega0=45, i=30, omega0=-40, theta0=theta0, t0_clock=t0_clock)

# Plotar a órbita
#plotter = Plotter(plot3d=True)
#plotter.plot_orbit(orbita, frame="bodycentric", points=True, velocities=True, positions=True)

#print(orbita)

maneuver = Maneuver.from_RTN(2, 0, 0, t_clock=1300, name="maneuver")

trajectory = Trajectory(orbita, name="Initial Orbit")
trajectory.add_maneuver(0, maneuver)
trajectory.add_trajectory_position(1, t_clock=5000, name="Final position")


# Plotar a órbita
plotter = Plotter(plot3d=True)
plotter.plot_trajectory(trajectory, frame="bodycentric",
                        points=True,
                        velocities=True,
                        positions=True,
                        orbits=True,
                        v_scale=5)

#print(trajectory)