from astronautica import Orbit, Body, Plotter, Trajectory
import numpy as np

earth = Body("earth")

t0_clock = 0
t1_clock = 4*24*3600

a = 26550
e = 0.74
i = np.degrees(np.arcsin(np.sqrt(2/5*2)))
Omega0 = 0
omega0 = 270
theta0 = 0

orbita = Orbit.from_elements(main_body=earth, a=a, e=e, i=i, Omega0=Omega0, omega0=omega0, theta0=theta0)


trajetoria = Trajectory(orbit0=orbita, t0_clock=t0_clock)
trajetoria.add_trajectory_position(0, t_clock=t1_clock, name="Final Position")

r_vec, v_vec = orbita.state_vectors_at_t_clock(t1_clock, frame="bodycentric")
print(r_vec)
print(v_vec)

plotter = Plotter(plot3d=True)
plotter.plot_trajectory(trajetoria, samples=1000, frame="rotating_bodycentric", orbits=False, velocities=True)

print(orbita)