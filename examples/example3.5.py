from orbitalmechanics import Orbit

"""
A geocentric trajectory has a perigee velocity of 15 km/s and a perigee altitude of 300 km. Find
(a) the radius and the time when the true anomaly is 100°;
(b) the position and speed 3 h later
"""

M_earth = 5.974e24 # [kg]
Earth_radius = 6378 # [km]
v = 15 # [km/s]
rp = 300 + Earth_radius # [km]
h = v * rp # [km^2/s]

orbita = Orbit(m1=M_earth, rp=rp, h=h, body1radius=Earth_radius)
orbita.add_position(0, "perigeu")

# (a)
theta0_a = 100
t1_a = orbita.t_from_theta(theta0_a)
r1_a = orbita.r_at_theta(theta0_a)
orbita.add_position(theta0_a, "100° após perigeu")

# (b)
t0_b = 3*60*60 + t1_a # [s]
theta1_b = orbita.theta_from_t(t0_b)
r1_b = orbita.r_at_theta(theta1_b)
v1_b = orbita.v_at_r(r1_b)
orbita.add_position(theta1_b, "3 horas após perigeu")

orbita.plot(plot_positions=True)

print(f"Tempo 100° após perigeu: {t1_a/(60*60):.4f} h")
print(f"Posição 100° após perigeu: {r1_a:.4f} km")
print(f"Posição 3 horas após perigeu: {r1_b:.4f} km")
print(f"Velocidade 3 horas após perigeu: {v1_b:.4f} km/s")




