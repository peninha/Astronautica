from orbitalmechanics import Orbit
import numpy as np

"""
At two points on a geocentric orbit, the altitude and true anomaly are z1 = 1545 km, θ1 = 126° and z2 = 852 km, θ2 = 58°,
respectively. Find (a) the eccentricity, (b) the altitude of perigee, (c) the semimajor axis, and (d) the period.
"""

M_earth = 5.9722e24 # [kg]
# Raio da Terra em km
R_terra = 6378.137

# Dados do problema
z1 = 1545  # km
theta1 = 126  # graus
z2 = 852   # km
theta2 = 58  # graus

# Convertendo altitudes para raios
r1 = R_terra + z1
r2 = R_terra + z2

# Criando uma órbita a partir dos dois pontos
orbita = Orbit.init_from_points(m1=M_earth, r1=r1, theta1=theta1, r2=r2, theta2=theta2, body1radius=R_terra)

# Obtendo os parâmetros orbitais
e = orbita.e
a = orbita.a
periodo = orbita.T/60  # convertendo para minutos
zp = orbita.rp - R_terra  # altitude do perigeu

# Imprimindo resultados
print(f"a) Excentricidade: {e:.4f}")
print(f"b) Altitude do perigeu: {zp:.2f} km")
print(f"c) Semi-eixo maior: {a:.2f} km")
print(f"d) Período: {periodo/60:.2f} horas")

# Plotando a órbita
orbita.plot(points=[(r1, theta1), (r2, theta2)])
