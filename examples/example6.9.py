from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
Spacecraft B and C are both in the geocentric elliptical orbit (1) shown in Fig. 6.20, from which it can be seen that the true
anomalies are θB ¼ 45° and θC ¼ 150°. At the instant shown, spacecraft B executes a delta-v maneuver, embarking upon a
trajectory (2), which will intercept and rendezvous with vehicle C in precisely one hour. Find the orbital parameters (e and
h) of the intercept trajectory and the total delta-v required for the chase maneuver.
"""
earth = Body("spherical_earth")

t0_clock = 0
rp = 8100
ra = 18900
orbit_from = Orbit.from_elements(earth, rp=rp, ra=ra, theta0=45, t0_clock=t0_clock)


orbit_of_target = Orbit.from_elements(earth, rp=rp, ra=ra, theta0=150, t0_clock=t0_clock)
theta_target = orbit_of_target.theta_at_t_clock(t0_clock+3600)

r_vec1, v_vec_A = orbit_from.state_vectors_at_t_clock(t0_clock)
r_vec2 = orbit_of_target.state_vectors_at_t_clock(t0_clock+3600)[0]

delta_t = 3600
orbit_to = Orbit.from_2_vectors_and_delta_time(earth, r_vec1, r_vec2, delta_t, t0_clock=t0_clock)

v_vec_B = orbit_to.state_vectors_at_t_clock(t0_clock-delta_t)[1]
delta_v_vec = v_vec_B - v_vec_A
maneuver = Maneuver.from_delta_v_vec(delta_v_vec, orbit_from, t0_clock, name="Maneuver", orbit_name="Intercept orbit", position_name="Burn position")

traj1 = Trajectory(orbit_from, orbit_name="Starting orbit", position_name="Burning position")
traj1.add_trajectory_position(0, t_clock=t0_clock-600, name="Starting position")
traj1.add_maneuver(0, maneuver, position_name="Burning position")
traj1.add_trajectory_position(1, t_clock=traj1.get_trajectory_position(0, position_index="last")['t_clock']+delta_t, name="Intercept position")

traj2 = Trajectory(orbit_of_target, orbit_name="Target orbit", position_name="Target at maneuver time")
traj2.add_trajectory_position(0, t_clock=t0_clock-600, name="Target starting position")
traj2.add_trajectory_position(0, t_clock=traj2.get_trajectory_position(0, position_index="last")['t_clock']+delta_t, name="Target intercept position")

plotter = Plotter(plot3d=False)
plotter.plot_trajectories([traj1, traj2],
                        frame="bodycentric",
                        points=True,
                        velocities=True,
                        positions=True,
                        orbits=False,
                        maneuvers=True,
                        time_step=60,
                        v_scale=1)
plotter.plot_orbit(orbit_of_target,
                    frame="bodycentric",
                    points=False,
                    velocities=False,
                    positions=False)

