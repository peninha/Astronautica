from astronautica import Orbit, Body, Plotter, Trajectory
import numpy as np

earth = Body("earth")

t0_clock = 0

rp = earth.radius_from_altitude(951)
e = 0.2421281
i = 44.8013
Omega0 = 170.9555
omega0 = 227.6939
theta0 = 0



orbita = Orbit.from_elements(main_body=earth, rp=rp, e=e, i=i, Omega0=Omega0, omega0=omega0, theta0=theta0)


t1_clock = 1*24*3600

trajetoria = Trajectory(orbit0=orbita)
trajetoria.add_trajectory_position(0, t_clock=t1_clock, name="Final Position")

r_vec, v_vec = orbita.state_vectors_at_t_clock(t1_clock, frame="bodycentric")
print(r_vec)
print(v_vec)

plotter = Plotter(plot3d=True)
plotter.plot_trajectory(trajetoria, frame="rotating_bodycentric", orbits=False, velocities=True, time_step=100)


print(orbita)