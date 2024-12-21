from orbitalmechanics import Orbit
import numpy as np

"""
Given the following state vector of a satellite in geocentric equatorial coordinates,
r = -3670I - 3870J + 4400K [km]
v = 4.7I - 7.4J + 1K [km/s]
find the state vector after 4 days (96 h) of coasting flight, assuming that there are no perturbations other than the influence of
the earth’s oblateness on Ω and ω.
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137
omega_earth = np.degrees(7.292115e-5) # [°/s]
J2 = 1.0826e-3
#J2 = 50*J2

t0_clock = 0
t1_clock = 96*3600

r_vec = np.array([-3670, -3870, 4400])
v_vec = np.array([4.7, -7.4, 1])

orbita = Orbit.init_from_state_vectors(r_vec, v_vec, m1=M_earth, m2=0, t0_clock=t0_clock, body1radius=R_terra)
orbita.add_oblateness_correction(J2, body1radius=R_terra)
orbita.add_main_body_rotation(omega_earth)

orbita.trajectory(t1_clock=t1_clock, n_points=100)

orbita.plot(frame="rotatingBodycentric",
            plot3d=True, 
            orbit=True,
            points=False,
            positions=True,
            velocities=True,
            trajectory=True)

r_vec, v_vec = orbita.state_vectors_at_t_clock(t1_clock, frame="bodycentric")
print(r_vec)
print(v_vec)