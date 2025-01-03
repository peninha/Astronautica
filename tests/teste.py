from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

earth = Body("squished_earth")

rp = 6800
ra = 13600
phase_change = -90
theta1_burn =45
n = 1

theta1_0 = 45
theta2_0 = theta1_0 + phase_change

orbita1 = Orbit.from_apsis(earth, rp, ra, theta0=theta1_0)
orbita2 = Orbit.from_apsis(earth, rp, ra, theta0=theta2_0)
traj1 = Trajectory(orbita1, name="Chaser")
traj2 = Trajectory(orbita2, name="Target")

t_burn = orbita1.t_clock_at_theta(theta1_burn)
theta2_burn = orbita2.theta_at_t_clock(t_burn)
phase_maneuver1, phase_maneuver2 = Maneuver.phase_maneuver(orbita1, theta2_burn-theta1_burn, theta_burn=theta1_burn, n=n)

traj1.add_maneuver(0, phase_maneuver1)
traj1.add_maneuver(1, phase_maneuver2)

t_maneuver1 = traj1.get_trajectory_position(0, position_index="last")['t_clock']
t_maneuver2 = traj1.get_trajectory_position(1, position_index="last")['t_clock']

traj2.add_trajectory_position(0, t_maneuver1, name="Target at Maneuver 1")
traj2.add_trajectory_position(0, t_maneuver2, name="Target at Maneuver 2")
#traj2.add_trajectory_position(0, t_maneuver2+1000, name="Target Final")
#traj1.add_trajectory_position(2, t_maneuver2+1000, name="Chaser Final")

plotter = Plotter(plot3d=False)
plotter.plot_trajectories([traj1, traj2],
                        frame="perifocal_t0",
                        points=True,
                        velocities=False,
                        positions=True,
                        orbits=False,
                        maneuvers=True,
                        time_step=60,
                        v_scale=1)

