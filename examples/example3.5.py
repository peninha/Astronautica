from astronautica import Orbit, Body, Plotter

"""
A geocentric trajectory has a perigee velocity of 15 km/s and a perigee altitude of 300 km. Find
(a) the radius and the time when the true anomaly is 100°;
(b) the position and speed 3 h later
"""

earth = Body(name="earth")
v = 15 # [km/s]
rp = 300 + earth.radius # [km]
h = v * rp # [km^2/s]

orbita = Orbit.from_elements(earth, rp=rp, h=h, theta0=0)

# (a)
theta0_a = 100
t1_a = orbita.t_clock_at_theta(theta0_a)
r1_a = orbita.r_at_theta(theta0_a)
orbita.add_orbital_position(theta=theta0_a, name="100° após perigeu")

# (b)
t0_b = 3*60*60 + t1_a # [s]
theta1_b = orbita.theta_at_t_clock(t0_b)
r1_b = orbita.r_at_theta(theta1_b)
v1_b = orbita.v_at_r(r1_b)
orbita.add_orbital_position(t_clock=t0_b, name="3 horas após perigeu")
    
plot = Plotter(frame="bodycentric", plot3d=True)
plot.plot_orbit(orbit=orbita)

print(f"Tempo 100° após perigeu: {t1_a/(60*60):.4f} h")
print(f"Posição 100° após perigeu: {r1_a:.4f} km")
print(f"Posição 3 horas após perigeu: {r1_b:.4f} km")
print(f"Velocidade 3 horas após perigeu: {v1_b:.4f} km/s")




