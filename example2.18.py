from orbitalmechanics import Three_body_restricted
import numpy as np

"""
For the spacecraft in the initial conditions (t = 0) are d = 200 km, ϕ = -90°, γ = 20°, and vbo = 10.9148 km/s.
Use the circular restricted three-body equations of motion, to determine the trajectory
and locate its position at t = 3.16689 days.    
"""

M_earth = 5.974e24 # [kg]
M_moon = 7.348e22 # [kg]
r12 = 3.844e5 # [km]
Earth_radius = 6378 # [km]
Moon_radius = 1737 # [km]
d = 200 # [km]
v = 10.9148 # [km/s]

three_body = Three_body_restricted(m1=M_earth, m2=M_moon, r12=r12, body1radius=Earth_radius, body2radius=Moon_radius)

r_vec_0 = np.array([-three_body.pi2*three_body.r12, -(d + Earth_radius), 0])
v_vec_0 = np.array([v*np.cos(np.radians(20)), -v*np.sin(np.radians(20)), 0])

sol = three_body.trajectory(r_vec_0, v_vec_0, t_span=(0, 3.16689*86400), t_eval=np.linspace(0, 3.16689*86400, 1000), plot=True)

# Obter posição final
r_final = sol.y[:3, -1]
print(f"\nPosição final (x,y,z) = ({r_final[0]:.2f}, {r_final[1]:.2f}, {r_final[2]:.2f}) km")
# Calcular a distância até a Lua
r_moon = np.array([(1-three_body.pi2)*three_body.r12, 0, 0])
dist_to_moon = np.linalg.norm(r_final - r_moon) - Moon_radius
print(f"\nAltitude em relação à Lua = {dist_to_moon:.2f} km")
