from astronautica import Orbit, Body, Plotter, Frame, Maneuver, Trajectory, Groundstations
import numpy as np

earth = Body(name="earth")

e = 0
rp = earth.radius_from_altitude(250)
theta0 = 7
t0_clock = 0
Omega0 = 0
i = -18.4
omega0 = 0

parking_orbit = Orbit.from_elements(earth, rp=rp, e=e, Omega0=Omega0, i=i, omega0=omega0, theta0=theta0, t0_clock=t0_clock)

traj = Trajectory(parking_orbit, orbit_name="Initial Orbit", position_name="Lançamento")
traj.add_trajectory_position(0, t_clock=400, name="Medição 1")
traj.add_trajectory_position(0, t_clock=5660, name="Medição 2")
traj.add_trajectory_position(0, t_clock=6060, name="Envio de manobra")


"""
# Plotar a órbita
plotter = Plotter(plot3d=True)
plotter.plot_trajectory(traj, frame="rotating_bodycentric",
                   points=True,
                   velocities=True,
                   positions=True,
                   time_step=60,
                   v_scale=1)
"""

stations = Groundstations(body=earth)
stations.add_groundstation(latitude=-2.316, longitude=-44.3676, name="Alcâncara")
stations.add_groundstation(latitude=-5.866, longitude=-35.383, name="Barreira do Inferno")
stations.add_groundstation(latitude=-23.218, longitude=-45.871, name="ITA")

frame = Frame.rotating_bodycentric(earth, Omega0=51, t0_clock=0)
Plotter.plot_groundtrack(traj, frame=frame, # type: ignore
                         earth_map=5,
                         groundstations=stations,
                         LOS=True,
                         time_step=60)

print(parking_orbit.T)

