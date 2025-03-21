from astronautica import Orbit, Body, Plotter, Trajectory
import numpy as np

earth = Body("earth")

t0_clock = 0
t1_clock = 16*24*3600

rp = earth.radius_from_altitude(300)
ra = earth.radius_from_altitude(20385)
i = np.degrees(np.arcsin(np.sqrt(2/5*2)))
Omega0 = 0
omega0 = 270
theta0 = 0

orbita = Orbit.from_elements(main_body=earth, rp=rp, ra=ra, i=i, Omega0=Omega0, omega0=omega0, theta0=theta0)


trajetoria = Trajectory(orbit0=orbita)
trajetoria.add_trajectory_position(0, t_clock=t1_clock, name="Final Position")

r_vec, v_vec = orbita.state_vectors_at_t_clock(t1_clock, frame="bodycentric")
print(r_vec)
print(v_vec)

plotter = Plotter(plot3d=True)
plotter.plot_trajectory(trajetoria, frame="rotating_bodycentric", orbits=False, velocities=True)

print(orbita)

print("Periodo: ", orbita.T / 3600, " horas")