from orbitalmechanics import Orbit
import numpy as np

"""
An earth satellite is in an orbit with a perigee altitude zp = 400 km and an apogee altitude za = 4000 km, as shown in
Fig. 2.21. Find each of the following quantities:
(a) eccentricity, e
(b) angular momentum, h
(c) perigee velocity, vp
(d) apogee velocity, va
(e) semimajor axis, a
(f) period of the orbit, T
(g) true anomaly-averaged radius rθ
(h) true anomaly when r ¼ rθ
(i) satellite speed when r ¼ rθ
(j) flight path angle γ when r ¼ rθ
(k) maximum flight path angle γmax and the true anomaly at which it occurs
"""

M_earth = 5.9722e24 # [kg]
# Raio da Terra em km
R_terra = 6378.137

# Altitudes do perigeu e apogeu
zp = 400  # km 
za = 4000  # km

# Distâncias do perigeu e apogeu ao centro da Terra
rp = R_terra + zp
ra = R_terra + za

# Criar órbita
orbita = Orbit(m1=M_earth, rp=rp, ra=ra, body1radius=R_terra)

# a) Excentricidade
print(f"\na) Excentricidade: {orbita.e:.4f}")

# b) Momento angular específico
print(f"b) Momento angular específico: {orbita.h:.2f} km²/s")

# c) Velocidade no perigeu
vp = orbita.v_at_r(orbita.rp)
print(f"c) Velocidade no perigeu: {vp:.3f} km/s")

# d) Velocidade no apogeu
va = orbita.v_at_r(orbita.ra)
print(f"d) Velocidade no apogeu: {va:.3f} km/s")

# e) Semi-eixo maior
print(f"e) Semi-eixo maior: {orbita.a:.2f} km")

# f) Período orbital
print(f"f) Período orbital: {orbita.T/60/60:.3f} horas")

# g) Raio médio
r_mean = orbita.r_avg_by_theta()
print(f"g) Raio médio: {r_mean:.2f} km")

# h) Anomalia verdadeira em raio médio
theta_r_mean = orbita.theta_at_r(r_mean)
print(f"h) Anomalia verdadeira em raio médio: {theta_r_mean[0]:.2f}° e {theta_r_mean[1]:.2f}°")

# i) Velocidade em raio médio
v_mean = orbita.v_at_r(r_mean)
print(f"i) Velocidade em raio médio: {v_mean:.3f} km/s")

# j) Ângulo de trajetória em raio médio
gamma_r_mean = orbita.gamma_at_theta(theta_r_mean[0])
print(f"j) Ângulo de trajetória em raio médio: {gamma_r_mean:.2f}°")

# k) Ângulo de trajetória máximo
theta_max = np.degrees(np.arccos(-orbita.e))
gamma_max = orbita.gamma_at_theta(theta_max)
print(f"k) Ângulo de trajetória máximo: {gamma_max:.2f}° em theta = {theta_max:.2f}°")

# Plotar a órbita
orbita.add_orbital_position(0, name="Perigeu")
orbita.add_orbital_position(180, name="Apogeu")
orbita.plot()