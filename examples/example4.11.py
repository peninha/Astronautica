from astronautica import Orbit, Body, Plotter, Trajectory
import numpy as np

"""
Given the following state vector of a satellite in geocentric equatorial coordinates,
r = -3670I - 3870J + 4400K [km]
v = 4.7I - 7.4J + 1K [km/s]
find the state vector after 4 days (96 h) of coasting flight, assuming that there are no perturbations other than the influence of
the earth’s oblateness on Ω and ω.
"""

earth = Body("earth")

t0_clock = 0
t1_clock = 96*3600

r_vec = np.array([-3670, -3870, 4400])
v_vec = np.array([4.7, -7.4, 1])

orbita = Orbit.from_state_vectors(earth, r_vec, v_vec, t0_clock=t0_clock)


trajetoria = Trajectory(orbit0=orbita, t0_clock=t0_clock)
trajetoria.add_trajectory_position(0, t_clock=t1_clock, name="Final Position")

r_vec, v_vec = orbita.state_vectors_at_t_clock(t1_clock, frame="bodycentric")
print(r_vec)
print(v_vec)

plotter = Plotter(plot3d=True)
plotter.plot_trajectory(trajetoria, samples=1000, frame="bodycentric", orbits=False, velocities=True)
