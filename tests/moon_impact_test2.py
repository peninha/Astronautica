from astronautica import Orbit, Body, Plotter, Frame, Maneuver, Trajectory
import numpy as np
import pandas as pd

earth = Body(name="spherical_earth")

r_vec = np.array([3.9019E+05, -7.6523E+04, -7.0725E+04])
v_vec = np.array([2.4873E-01, 8.7246E-01, 3.4007E-01])
t0_clock = 0

moon_orbit = Orbit.from_state_vectors(earth, r_vec, v_vec, t0_clock)

"""
df = pd.read_csv("../data/moon_vectors_01.csv")
for i in range(24):
    t = i * 3600
    theta = moon_orbit.theta_at_t_clock(t)
    r_vec, v_vec = moon_orbit.state_vectors_at_t_clock(t, frame="bodycentric")
    
    # Carrega os dados do arquivo CSV
    x_csv = df.iloc[i, 2]  # Coluna C (índice 2)
    y_csv = df.iloc[i, 3]  # Coluna D (índice 3) 
    z_csv = df.iloc[i, 4]  # Coluna E (índice 4)
    vx_csv = df.iloc[i, 5]  # Coluna F (índice 5)
    vy_csv = df.iloc[i, 6]  # Coluna G (índice 6)
    vz_csv = df.iloc[i, 7]  # Coluna H (índice 7)
    
    r_vec_csv = np.array([x_csv, y_csv, z_csv])
    v_vec_csv = np.array([vx_csv, vy_csv, vz_csv])
    
    print(f"\n#########\nt = {i} horas")
    print(f"r_vec calculado = {r_vec}")
    print(f"r_vec do arquivo = {r_vec_csv}")
    print(f"diferença r_vec = {r_vec - r_vec_csv}")
    print(f"v_vec calculado = {v_vec}")
    print(f"v_vec do arquivo = {v_vec_csv}")
    print(f"diferença v_vec = {v_vec - v_vec_csv}")

"""


e1 = 0
rp1 = earth.radius_from_altitude(250)
i1 = moon_orbit.i -1
Omega01 = moon_orbit.Omega0
omega01 = 0
theta01 = 0
t0_clock1 = 0

parking_orbit = Orbit.from_elements(earth, e=e1, rp=rp1, Omega0=Omega01, i=i1, omega0=omega01, theta0=theta01, t0_clock=t0_clock1)

theta_rel = 138
from_t_clock = 300
delta_v_target = 3.2
delta_t_min = 1.5*24*3600
delta_t_max = 3.5*24*3600
delta_t_guess = 2.5*24*3600


"""
t_clock_burn = Maneuver.find_t_clock_for_relative_theta(parking_orbit, moon_orbit, theta_rel, from_t_clock=from_t_clock, max_iter=1000, tol=1e-6)
RPN = np.array([0, delta_v_target, 0])
burn = Maneuver.from_RPN(RPN, t_clock_burn, name="Translunar injection")

t_clock_impact = t_clock_burn + delta_t_guess

traj1 = Trajectory(parking_orbit, orbit_name="Parking orbit", position_name="Probe initial position")
traj1.add_maneuver(0, burn, name="Translunar injection")
traj1.add_trajectory_position(1, t_clock=t_clock_impact, name="Probe impact position")

traj2 = Trajectory(moon_orbit, orbit_name="Moon orbit", position_name="Moon initial position")
traj2.add_trajectory_position(0, t_clock=t_clock_burn, name="Moon at burn")
traj2.add_trajectory_position(0, t_clock=t_clock_impact, name="Moon at impact")

plotter = Plotter(plot3d=True)
plotter.plot_trajectories([traj1, traj2], frame="perifocal",
                           time_step=600,
                           orbits=False)

v_impact_vec = traj1.orbits[1]['orbit'].state_vectors_at_t_clock(t_clock_impact, frame="perifocal")[1]
v_impact = np.linalg.norm(v_impact_vec)
print("v_impact_vec: ", v_impact_vec)
print("v_impact: ", v_impact)
"""

result = Maneuver.impact_maneuver(parking_orbit, moon_orbit, theta_rel, delta_v_target, delta_t_guess, delta_t_min, delta_t_max, from_t_clock=from_t_clock, tol=1e-6, max_iter=1000)

burn = result["maneuver"]
delta_t = result["delta_t"]
delta_v = result["delta_v"]
t_clock_burn = result["t_clock_burn"]
t_clock_impact = t_clock_burn + delta_t

traj1 = Trajectory(parking_orbit, orbit_name="Parking orbit", position_name="Probe initial position")
traj1.add_maneuver(0, burn, name="Translunar injection")
traj1.add_trajectory_position(1, t_clock=t_clock_impact, name="Probe impact position")

traj2 = Trajectory(moon_orbit, orbit_name="Moon orbit", position_name="Moon initial position")
traj2.add_trajectory_position(0, t_clock=t_clock_burn, name="Moon at burn")
traj2.add_trajectory_position(0, t_clock=t_clock_impact, name="Moon at impact")

plotter = Plotter(plot3d=True)
plotter.plot_trajectories([traj1, traj2], frame="bodycentric",
                           time_step=600,
                           orbits=False)


v_impact_vec_ship = traj1.orbits[1]['orbit'].state_vectors_at_t_clock(t_clock_impact, frame="bodycentric")[1]
v_impact_vec_moon = traj2.orbits[0]['orbit'].state_vectors_at_t_clock(t_clock_impact, frame="bodycentric")[1]
v_impact_vec = v_impact_vec_ship - v_impact_vec_moon
v_impact = np.linalg.norm(v_impact_vec)

print("delta_t: ", delta_t)
print("delta_v: ", delta_v)
print("t_clock_burn: ", f"{int(t_clock_burn//3600)}h {int((t_clock_burn%3600)//60)}m {int(t_clock_burn%60)}s")
print("t_clock_impact: ", f"{int(t_clock_impact//3600)}h {int((t_clock_impact%3600)//60)}m {int(t_clock_impact%60)}s")
print("v_impact_vec_ship: ", v_impact_vec_ship)
print("v_impact_vec_moon: ", v_impact_vec_moon)
print("v_impact_vec: ", v_impact_vec)
print("v_impact: ", v_impact)
print("RTN:", burn.RTN)
