from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
A geocentric satellite in orbit 1 of Fig. 6.15 executes a delta-v maneuver at A,
which places it on orbit 2, for reentry at D.
Calculate Δv at A and its direction relative to the local horizon
"""
earth = Body("spherical_earth")

rp = 10000
ra = 20000
theta0 = 150
t0_clock = 0
orbit_from = Orbit.from_elements(earth, rp=rp, ra=ra, theta0=theta0, t0_clock=t0_clock)


r1 = orbit_from.r_at_theta(theta0)
theta1 = theta0
r2 = earth.radius_from_altitude(0)
theta2 = 0
orbit_to = Orbit.from_2_positions(earth, r1, theta1, r2, theta2, t0_clock=0)

delta_v_vec = orbit_to.state_vectors_at_theta(theta1)[1] - orbit_from.state_vectors_at_theta(theta0)[1]
maneuver = Maneuver.from_delta_v_vec(delta_v_vec, orbit_from, t0_clock)

if maneuver is None:
    print("No maneuver found")
else:
    if maneuver.RTN is not None:
        angle_to_local_horizon = np.degrees(np.arctan2(maneuver.RTN[0], maneuver.RTN[1]))
        print(f"Delta-v: {np.linalg.norm(delta_v_vec)} km/s")
        print(f"Angle to local horizon: {angle_to_local_horizon} deg")
    else:
        print("Maneuver is not in RTN mode")

    traj1 = Trajectory(orbit_from, orbit_name="Starting orbit", position_name="Burn position")
    traj1.add_trajectory_position(0, t_clock=t0_clock-1000, name="Starting position")
    traj1.add_maneuver(0, maneuver, position_name="Maneuver position")
    traj1.add_trajectory_position(1, t_clock=-orbit_to.t_clock_at_theta(-360+theta1), name="Ground position")

    plotter = Plotter(plot3d=False)
    plotter.plot_trajectory(traj1,
                            frame="bodycentric",
                            points=True,
                            velocities=True,
                            positions=True,
                            orbits=True,
                            maneuvers=True,
                            time_step=300,
                            v_scale=5)

