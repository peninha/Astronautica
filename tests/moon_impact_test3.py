from astronautica import Orbit, Body, Plotter, Frame, Maneuver, Trajectory
import numpy as np
import pandas as pd

earth = Body(name="earth")


rp = earth.radius_from_altitude(245.52336)
ra = earth.radius_from_altitude(252.0211)
i = 19.90664
Omega0 = 349.79269
omega0 = 312.58068
theta0 = 251.63458
t0_clock = 0

parking_orbit = Orbit.from_elements(earth, rp=rp, ra=ra, Omega0=Omega0, i=i, omega0=omega0, theta0=theta0, t0_clock=t0_clock)

theta_rel = 138
from_t_clock = 300
delta_v_target = 3.224
delta_t_min = 1.0*24*3600
delta_t_guess = 2.0*24*3600
delta_t_max = 3.0*24*3600




rp = earth.radius_from_altitude(352601.66175)
ra = earth.radius_from_altitude(406060.52059)
print(rp, ra)

i = 19.91161
Omega0 = 349.79818
omega0 = 180.84834
theta0 = 293.10453
t0_clock = 0

moon_orbit = Orbit.from_elements(earth, rp=rp, ra=ra, Omega0=Omega0, i=i, omega0=omega0, theta0=theta0, t0_clock=t0_clock)

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
                           orbits=True)


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
print("TNB:", parking_orbit.convert_RTN_to_TNB(burn.RTN, t_clock=t_clock_burn))
print("Phase angle of t0_clock:", Maneuver.phase_angle_at_t_clock(parking_orbit, moon_orbit, t0_clock))
