from astronautica import Orbit, Body, Plotter, Frame, Maneuver, Trajectory
import numpy as np
import pandas as pd

earth = Body(name="earth")

e = 6.253105788238640E-02
rp = 3.570157160195780E+05
i = 2.325318631898455E+01
Omega0 = 1.335674839291536E+01
omega0 = 1.675204123831447E+02
theta0 = 1.661472632534166E+02
t0_clock = 0

moon_orbit = Orbit.from_elements(earth, e=e, rp=rp, Omega0=Omega0, i=i, omega0=omega0, theta0=theta0, t0_clock=t0_clock)

df = pd.read_csv("../data/moon_01.csv")
for i in range(24):
    t = i * 3600
    theta = moon_orbit.theta_at_t_clock(t)
    r_vec, v_vec = moon_orbit.state_vectors_at_t_clock(t, frame="bodycentric")
    
    # Carrega os dados do arquivo CSV
    theta_csv = df.iloc[i, 10]  # Coluna 11 (índice 10)
    
    print(f"\n#########\nt = {i} horas")
    print(f"theta calculado = {theta}")
    print(f"theta do arquivo = {theta_csv}")
    print(f"diferença = {theta - theta_csv}")
    print(f"r_vec = {r_vec}")
    print(f"v_vec = {v_vec}")



"""
e1 = 0
rp1 = earth.radius_from_altitude(200)
i1 = i
Omega01 = 0
omega01 = 0
theta01 = 0
t0_clock1 = 0

orbit_from = Orbit.from_elements(earth, e=e1, rp=rp1, Omega0=Omega01, i=i1, omega0=omega01, theta0=theta01, t0_clock=t0_clock1)

plotter = Plotter(plot3d=True)
plotter.plot_orbits([moon_orbit, orbit_from], frame="bodycentric", h_arrow=False, e_arrow=True, n_arrow=False)
"""