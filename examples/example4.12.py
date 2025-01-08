from astronautica import Orbit, Body, Frame, Plotter, Trajectory
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

earth = Body("earth")

rp = 6700 # [km]
ra = 10000 # [km]
Omega0 = 270 # [°]
i0 = 60 # [°]
omega0 = 45 # [°]
omega_earth = np.degrees(7.292115e-5) # [°/s]

theta0 = 230 # [°]
t0_clock = 0 # [s]
t1_clock = 45*60*1 # [s]

orbita = Orbit.from_elements(main_body=earth, rp=rp, ra=ra, theta0=theta0, Omega0=Omega0, i=i0, omega0=omega0)
print(orbita)

theta1 = orbita.theta_at_t_clock(t1_clock)
print("theta0: ", theta0)
print("theta1: ", theta1)

r_vec_bc_1, v_vec_bc_1 = orbita.state_vectors_at_t_clock(t1_clock, frame="bodycentric")
print("r_vec_bc_1: ", r_vec_bc_1)

r_vec_rbc_1, v_vec_rbc_1 = orbita.state_vectors_at_t_clock(t1_clock, frame="rotating_bodycentric")
print("r_vec_rbc_1: ", r_vec_rbc_1)

ra_1, dec_1 = Frame.convert_cartesian_to_ra_dec(r_vec_rbc_1)
print("ra_1: ", ra_1)
print("dec_1: ", dec_1)

trajetoria = Trajectory(orbita)
trajetoria.add_trajectory_position(0, t1_clock, theta1, name="Position 1")

plotter = Plotter(plot3d=True)
plotter.plot_trajectory(trajetoria,
                        frame="rotating_bodycentric",
                        time_step=60,
                        orbits=True,
                        points=True,
                        velocities=True,    
                        positions=True,
                        groundtrack=False)

plotter.plot_groundtrack(trajetoria, time_step=60)