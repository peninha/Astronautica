from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
It is desired to shift the longitude of a GEO satellite 12° westward in three revolutions of its phasing orbit.
Calculate the delta-v requirement.
"""
earth = Body("earth")

e = 0
rp = earth.geosynchronous_radius()
phase_change = -12
n = 3

orbita1 = Orbit.from_elements(earth, e=e, rp=rp, theta0=0)

delta_v1, delta_v2, _ = Maneuver.delta_v_for_phase_change(orbita1, phase_change, theta_burn=0, n=n)

phase_maneuver1, phase_maneuver2 = Maneuver.phase_maneuver(orbita1, phase_change, theta_burn=0, n=n)
traj1 = Trajectory(orbita1, orbit_name="Starting orbit", position_name="Chaser starting position")
traj1.add_maneuver(0, phase_maneuver1)
traj1.add_maneuver(1, phase_maneuver2)

print(f"Delta-v 1: {delta_v1} km/s")
print(f"Delta-v 2: {delta_v2} km/s")
print(f"Total delta-v: {np.abs(delta_v1) + np.abs(delta_v2)} km/s")

orbita2 = Orbit.from_elements(earth, e=e, rp=rp, theta0=phase_change)
traj2 = Trajectory(orbita2, orbit_name="Target orbit", position_name="Target starting position")
traj2.add_trajectory_position(0, traj1.get_trajectory_position(2, position_index="last")['t_clock'], name="Target end position")

plotter = Plotter(plot3d=False)
plotter.plot_trajectories([traj1, traj2],
                        frame="bodycentric",
                        points=True,
                        velocities=False,
                        positions=True,
                        orbits=False,
                        maneuvers=True,
                        time_step=600,
                        v_scale=30)

