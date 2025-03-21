from astronautica import Orbit, Body, Plotter, Frame, Maneuver, Trajectory
import numpy as np

earth = Body(name="earth")

rp = 7000
ra = 35000
theta0 = -250
t0_clock = 0

orbita = Orbit.from_elements(earth, rp=rp, ra=ra, Omega0=45, i=30, omega0=-40, theta0=theta0, t0_clock=t0_clock)

trajectory = Trajectory(orbita, orbit_name="Initial Orbit", position_name="Initial position")

maneuver = Maneuver.from_TNB([1, 0, 0], t_clock=4300, name="maneuver1")

trajectory.add_maneuver(0, maneuver)

maneuver2 = Maneuver.from_TNB([0, 1, 0], t_clock=50000, name="maneuver2")

trajectory.add_maneuver(1, maneuver2)

maneuver3 = Maneuver.from_TNB([0, 0, 1], t_clock=100000, name="maneuver3")

trajectory.add_maneuver(2, maneuver3)

maneuver4 = Maneuver.from_TNB([-1, 0, 0], t_clock=140000, name="maneuver4")

trajectory.add_maneuver(3, maneuver4)

trajectory.add_trajectory_position(4, t_clock=170000, name="Final position")

# Plotar a órbita
plotter = Plotter(plot3d=True)
plotter.plot_trajectory(trajectory, frame="bodycentric", time_step=600,
                        points=True,
                        velocities=True,
                        positions=True,
                        orbits=False,
                        maneuvers=True,
                        v_scale=5)

#print(trajectory)