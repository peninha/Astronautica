from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

earth = Body("squished_earth")

rp = 6800
ra = 13600
phase_change = -90
theta1_burn =45
n = 1

theta0 = 45
orbita0 = Orbit.from_apsis(earth, rp, ra, theta0=theta0)
traj0 = Trajectory(orbita0, name="Orbit 0")
traj0.add_trajectory_position(0, 1000)


v1 = np.float64(8.51500622059919)
a1 = np.float64(11027.619325681027)
h1 = np.float64(61370.35204448525)
e1 = np.float64(0.37836897972733635)

#gamma1 = orbita0.gamma_at_theta(theta1_burn)
#r = orbita0.r_at_theta(theta1_burn)
#v1 = orbita0.v_at_theta(theta1_burn)

r_vec, v_vec = orbita0.state_vectors_at_t_clock(0, frame="bodycentric")
v = np.linalg.norm(v_vec)
v1_vec = v_vec*v1/v

orbita1 = Orbit.from_state_vectors(earth, r_vec, v1_vec, t0_clock=0)
traj1 = Trajectory(orbita1, name="Orbit 1")
traj1.add_trajectory_position(0, 1000)



plotter = Plotter(plot3d=False)
plotter.plot_trajectories([traj0, traj1],
                        frame="perifocal_t0",
                        points=True,
                        velocities=False,
                        positions=True,
                        orbits=False,
                        maneuvers=True,
                        time_step=60,
                        v_scale=1)

