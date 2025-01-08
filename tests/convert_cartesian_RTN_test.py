from astronautica import Orbit, Body, Plotter, Frame, Maneuver, Trajectory
import numpy as np

earth = Body(name="spherical_earth")

rp = 10000
ra = 10000
Omega0 = 90
i = 60
omega0 = 0
theta0 = 90

orbit = Orbit.from_elements(earth, rp=rp, ra=ra, Omega0=Omega0, i=i, omega0=omega0, theta0=theta0, t0_clock=0)

#v_vec = orbit.state_vectors_at_theta(theta0, frame="bodycentric")[1]

#RTN = orbit.convert_cartesian_to_RTN(v_vec, 0, frame="bodycentric")

RTN = np.array([0, 0, 1])

maneuver = Maneuver.from_RTN(RTN, 0)

traj = Trajectory(orbit)

traj.add_maneuver(0, maneuver)

v_vec_peri = orbit.convert_RTN_to_cartesian(RTN, 0, frame="perifocal")
v_vec_bc = orbit.convert_RTN_to_cartesian(RTN, 0, frame="bodycentric")

print("RTN: ", RTN)
print("v_vec_peri: ", np.round(v_vec_peri, 2))
print("v_vec_bc: ", np.round(v_vec_bc, 2))

RTN1 = orbit.convert_cartesian_to_RTN(v_vec_peri, 0, frame="perifocal")
RTN2 = orbit.convert_cartesian_to_RTN(v_vec_bc, 0, frame="bodycentric")

print("RTN1: ",  np.round(RTN1, 2))
print("RTN2: ",  np.round(RTN2, 2))

#RTN1 = orbit_from.convert_cartesian_to_RTN(delta_v1_vec, t_clock_burn, frame="bodycentric")
plotter = Plotter(plot3d=True)
#plotter.plot_orbit(orbit, frame="bodycentric", points=True, velocities=True, positions=True)
plotter.plot_trajectory(traj, frame="bodycentric",
                        velocities=False,
                        positions=True,
                        maneuvers=True,
                        v_scale=10)

