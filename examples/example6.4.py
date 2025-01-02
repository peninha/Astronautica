from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
Spacecraft at A and B are in the same orbit (1). At the instant shown in Fig. 6.11 the chaser vehicle
at A executes a phasing maneuver so as to catch the target spacecraft back at A after just one revolution
of the chaser’s phasing orbit (2). What is the required total delta-v?
"""
earth = Body("earth")

rp = 6800
ra = 13600
phase_change = 90

orbita1 = Orbit.from_apsis(earth, rp, ra, theta0=0)

delta_v1, delta_v2, _ = Maneuver.delta_v_for_phase_change(orbita1, phase_change, theta_burn=0, n=1)

phase_maneuver1, phase_maneuver2 = Maneuver.phase_maneuver(orbita1, phase_change, theta_burn=0, n=1)
traj1 = Trajectory(orbita1)
traj1.add_maneuver(0, phase_maneuver1)
traj1.add_maneuver(1, phase_maneuver2)

print(f"Delta-v 1: {delta_v1} km/s")
print(f"Delta-v 2: {delta_v2} km/s")
print(f"Total delta-v: {np.abs(delta_v1) + np.abs(delta_v2)} km/s")

orbita2 = Orbit.from_apsis(earth, rp, ra, theta0=phase_change)
traj2 = Trajectory(orbita2)
traj2.add_trajectory_position(0, traj1.get_trajectory_position(2, position_index="last")['t_clock'], name="Position 2")


plotter = Plotter(plot3d=True)
plotter.plot_trajectories([traj1, traj2],
                        frame="bodycentric",
                        points=True,
                        velocities=False,
                        positions=True,
                        orbits=False,
                        maneuvers=True,
                        time_step=300,
                        v_scale=30)

