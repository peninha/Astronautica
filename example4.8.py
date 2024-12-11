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
theta0 = 230 # [°]
Omega0 = 270 # [°]
i0 = 60 # [°]
omega0 = 45 # [°]
omega_earth = np.degrees((2*np.pi + 2*np.pi/365.25)/(24*60*60)) # [°/s]

theta1 = 300 # [°]

orbita = Orbit(m1=M_earth, m2=0, rp=rp, ra=ra, theta=theta0, Omega=Omega0, i=i0, omega=omega0, body1radius=R_terra)
orbita.omega_body = omega_earth

r_vec_rbc_0, v_vec_rbc_0 = orbita.state_vectors_at_theta(theta0, frame="rotatingBodycentric")
r_vec_rbc_1, v_vec_rbc_1 = orbita.state_vectors_at_theta(theta1, frame="rotatingBodycentric")

print(r_vec_rbc_0)
print(r_vec_rbc_1)

ra_0, dec_0 = orbita.convert_cartesian_to_ra_dec(r_vec_rbc_0)
ra_1, dec_1 = orbita.convert_cartesian_to_ra_dec(r_vec_rbc_1)

print(ra_0, dec_0)
print(ra_1, dec_1)

#r_vec, v_vec = orbita.state_vectors_at_t(45*60, frame="rotatingBodycentric")


#orbita.plot(frame="bodycentric", points=True, velocities=True, positions=True, trajectory=False, plot3d=True)

