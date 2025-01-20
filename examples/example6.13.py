from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
A spacecraft is in a 500 km by 10,000 km altitude geocentric orbit that intersects the equatorial
plane at a true anomaly of 120° (see Fig. 6.33). If the orbit’s inclination to the equatorial plane is 15°,
what is the minimum velocity increment required to make this an equatorial orbit?
"""
earth = Body("spherical_earth")

rp = earth.radius_from_altitude(500)
ra = earth.radius_from_altitude(10000)

i1 = 15
Omega01 = 120 - 180
omega01 = 180 - 120

i2 = 0
Omega02 = 0
omega02 = 0

orbit_from = Orbit.from_elements(earth, rp=rp, ra=ra, i=i1, Omega0=Omega01, omega0=omega01, theta0=120, t0_clock=0, name="Orbit from")
orbit_to = Orbit.from_elements(earth, rp=rp, ra=ra, i=i2, Omega0=0, omega0=0, theta0=120, t0_clock=0, name="Orbit to")

maneuver, theta_burn, delta_v = Maneuver.orbit_change_maneuver(orbit_from, orbit_to, theta_burn=120)

traj1 = Trajectory(orbit_from)
traj1.add_trajectory_position(0, -orbit_from.T/2)
traj1.add_maneuver(0, maneuver)
traj1.add_trajectory_position(1, orbit_to.T/2)

print(f"Delta-v required: {delta_v:.5f} km/s")

plotter = Plotter(plot3d=True)
plotter.plot_trajectory(traj1, orbits=False)