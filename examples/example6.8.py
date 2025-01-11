from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
An earth satellite in orbit 1 of Fig. 6.19 undergoes the indicated delta-v maneuver at its perigee.
Determine the rotation η of its apse line as well as the new perigee and apogee.
"""
earth = Body("spherical_earth")

t0_clock = 0
rp = 7000
ra = 17000
orbit_from = Orbit.from_elements(earth, rp=rp, ra=ra, theta0=0, t0_clock=t0_clock)

delta_v = 2 # km/s
angle = 60 # deg
maneuver = Maneuver.from_delta_v_and_angle(delta_v, angle, t0_clock)


traj1 = Trajectory(orbit_from, orbit_name="Starting orbit", position_name="Burning position")
traj1.add_trajectory_position(0, t_clock=t0_clock-1000, name="Starting position")
traj1.add_maneuver(0, maneuver, position_name="Burning position")
traj1.add_trajectory_position(1, t_clock=traj1.get_trajectory_position(0, position_index="last")['t_clock']+10000, name="Final position")

eta = Orbit.convert_cartesian_to_polar(traj1.orbits[1]['orbit'].e_vec_bc)[1] - 360
print(f"Rotation η: {eta} deg")
print(f"Perigee: {traj1.orbits[1]['orbit'].rp} km")
print(f"Apogee: {traj1.orbits[1]['orbit'].ra} km")

plotter = Plotter(plot3d=False)
plotter.plot_trajectory(traj1,
                        frame="bodycentric",
                        points=True,
                        velocities=True,
                        positions=True,
                        orbits=True,
                        maneuvers=True,
                        time_step=300,
                        v_scale=3)
