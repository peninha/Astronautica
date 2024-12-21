from astronautica import Orbit, Body, Plotter
import numpy as np

"""
At two points on a geocentric orbit, the altitude and true anomaly are z1 = 1545 km, θ1 = 126° and z2 = 852 km, θ2 = 58°,
respectively. Find (a) the eccentricity, (b) the altitude of perigee, (c) the semimajor axis, and (d) the period.
"""

earth = Body(name="earth")

# Dados do problema
z1 = 1545  # km
theta1 = 126  # graus
z2 = 852   # km
theta2 = 58  # graus

# Convertendo altitudes para raios
r1 = earth.radius + z1
r2 = earth.radius + z2

# Criando uma órbita a partir dos dois pontos
orbita = Orbit.from_2_positions(earth, r1=r1, theta1=theta1, r2=r2, theta2=theta2)

# Obtendo os parâmetros orbitais
e = orbita.e
a = orbita.a
periodo = orbita.T/60  # convertendo para minutos
zp = orbita.rp - earth.radius  # altitude do perigeu

# Imprimindo resultados
print(f"a) Excentricidade: {e:.4f}")
print(f"b) Altitude do perigeu: {zp:.2f} km")
print(f"c) Semi-eixo maior: {a:.2f} km")
print(f"d) Período: {periodo/60:.2f} horas")

# Plotando a órbita
Plotter(frame="bodycentric", plot3d=True).plot_orbit(orbita)
