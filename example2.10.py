from orbitalmechanics import Orbit
import numpy as np

"""
At a given point of a spacecraft’s geocentric trajectory, the radius is 14,600 km, the speed is 8.6 km/s, and the flight path
angle is 50°. Show that the path is a hyperbola and calculate the following:
(a) angular momentum
(b) eccentricity
(c) true anomaly
(d) radius of the perigee
(e) semimajor axis
(f) C3
(g) turn angle
(h) aiming radius
"""
M_earth = 5.9722e24 # [kg]
# Raio da Terra em km
R_terra = 6378.137
# Dados do problema
r = 14600  # km
v = 8.6    # km/s
gamma = 50 # graus

orbita = Orbit.init_from_r_v_gamma(mu=M_earth, r=r, v=v, gamma=gamma, body1radius=R_terra)

# Calcular velocidade de escape na posição r
v_esc = orbita.v_esc(r)

# Verificar se é hiperbólica
print(orbita.type())

# a) Momento angular específico
print(f"\na) Momento angular específico: {orbita.h:.2f} km²/s")

# b) Excentricidade
print(f"b) Excentricidade: {orbita.e:.4f}")

# c) Anomalia verdadeira
theta = orbita.theta_at_r(r)
print(f"c) Anomalia verdadeira: {theta[0]:.2f}°")

# d) Raio do perigeu
print(f"d) Raio do perigeu: {orbita.rp:.2f} km")
print(f"   Altitude do perigeu: {orbita.rp - R_terra:.2f} km")

# e) Semi-eixo maior
print(f"e) Semi-eixo maior: {orbita.a:.2f} km")

# f) Energia característica (C3)
C3 = orbita.C3()
print(f"f) C3: {C3:.2f} km²/s²")

# g) Ângulo de virada
delta = orbita.turn_angle()
print(f"g) Ângulo de virada: {delta:.2f}°")

# h) Raio de mira
b = orbita.aiming_radius()
print(f"h) Raio de mira: {b:.2f} km")

# Plotar a órbita
orbita.plot(points=[(r, theta[0])])
