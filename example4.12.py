from orbitalmechanics import Orbit
import numpy as np

"""
An earth satellite has the following orbital parameters:
rp = 6700 km Perigee
ra = 10,000 km Apogee
θ0 = 230° True anomaly
Ω0 = 270° Right ascension of the ascending node
i0 = 60° Inclination
ω0 = 45° Argument of perigee
Calculate the right ascension (longitude east of x') and declination (latitude)
relative to the rotating earth 45 min later.
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137

rp = 6700 # [km]
ra = 10000 # [km]
Omega0 = 270 # [°]
i0 = 60 # [°]
omega0 = 45 # [°]
omega_earth = np.degrees((2*np.pi + 2*np.pi/365.25)/(24*60*60)) # [°/s]

theta0 = 230 # [°]
t0_clock = 0 # [s]
t1_clock = 45*60 # [s]

orbita = Orbit(m1=M_earth, m2=0, rp=rp, ra=ra, theta0=theta0, Omega=Omega0, i=i0, omega=omega0, body1radius=R_terra)
orbita.omega_body = omega_earth

#r_vec_rbc_0, v_vec_rbc_0 = orbita.state_vectors_at_theta(theta0, frame="rotatingBodycentric")
#r_vec_rbc_1, v_vec_rbc_1 = orbita.state_vectors_at_theta(theta1, frame="rotatingBodycentric")

#print(r_vec_rbc_0)
#print(r_vec_rbc_1)

#ra_0, dec_0 = orbita.convert_cartesian_to_ra_dec(r_vec_rbc_0)
#ra_1, dec_1 = orbita.convert_cartesian_to_ra_dec(r_vec_rbc_1)

#print(ra_0, dec_0)
#print(ra_1, dec_1)

delta_t = t1_clock - t0_clock
t0_orbit = orbita.t_orbit_at_theta(theta0)
t1_orbit = t0_orbit + delta_t
theta1 = orbita.theta_at_t_clock(t1_clock)
print("theta0: ", theta0)
print("theta1: ", theta1)

#r_vec_bc_0, v_vec_bc_0 = orbita.state_vectors_at_t_clock(t0_clock, frame="bodycentric")
r_vec_bc_1, v_vec_bc_1 = orbita.state_vectors_at_t_clock(t1_clock, frame="bodycentric")

#print("r_vec_bc_0: ", r_vec_bc_0)
print("r_vec_bc_1: ", r_vec_bc_1)

#r_vec_rbc_0, v_vec_rbc_0 = orbita.state_vectors_at_t_clock(t0_clock, frame="rotatingBodycentric")
r_vec_rbc_1, v_vec_rbc_1 = orbita.state_vectors_at_t_clock(t1_clock, frame="rotatingBodycentric")

#print("r_vec_rbc_0: ", r_vec_rbc_0)
print("r_vec_rbc_1: ", r_vec_rbc_1)

#ra_0, dec_0 = orbita.convert_cartesian_to_ra_dec(r_vec_rbc_0)
ra_1, dec_1 = orbita.convert_cartesian_to_ra_dec(r_vec_rbc_1)

#print("ra_0: ", ra_0)
#print("dec_0: ", dec_0)
print("ra_1: ", ra_1)
print("dec_1: ", dec_1)

orbita.add_zero_state_from_theta(theta0, t0_clock=t0_clock)

orbita.plot(frame="rotatingBodycentric", points=True, velocities=True, positions=True, trajectory=True, plot3d=True)

