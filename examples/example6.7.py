from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
An earth satellite is in an 8000 km by 16,000 km radius orbit (orbit 1 of Fig. 6.18).
Calculate the delta-v and the true anomaly θ1 required to obtain a
7000 km by 21,000 km radius orbit (orbit 2) whose apse line is rotated 25° counterclockwise.
Indicate the orientation ϕ of Δv to the local horizon.
"""
earth = Body("spherical_earth")

t0_clock = 0
rp1 = 8000
ra1 = 16000
orbit_from = Orbit.from_elements(earth, rp=rp1, ra=ra1, theta0=0, t0_clock=t0_clock)

rp2 = 7000
ra2 = 21000
omega0 = 25
orbit_to = Orbit.from_elements(earth, rp=rp2, ra=ra2, omega0=omega0, t0_clock=t0_clock)

maneuver, theta_burn, delta_v = Maneuver.orbit_change_maneuver(orbit_from, orbit_to)

print(f"Delta-v: {delta_v} km/s")    
print(f"Theta burn: {theta_burn} deg")

if maneuver is None:
    print("No maneuver found")
else:
    if maneuver.RTN is not None:
        angle_to_local_horizon = np.degrees(np.arctan2(maneuver.RTN[0], maneuver.RTN[1]))
        print(f"Angle to local horizon: {angle_to_local_horizon} deg")
    else:
        print("Maneuver is not in RTN mode")

    traj1 = Trajectory(orbit_from, orbit_name="Starting orbit", position_name="Starting position")
    traj1.add_maneuver(0, maneuver, position_name="Maneuver position")
    traj1.add_trajectory_position(1, t_clock=traj1.get_trajectory_position(0, position_index="last")['t_clock']+10000, name="Final position")

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

    #plotter.plot_orbits([orbit'_from, orbit_to], frame="bodycentric", points=True, velocities=True, positions=True)