from orbitalmechanics import Orbit
import numpy as np

"""
Calculate the altitude zGEO and speed vGEO of a geostationary earth satellite.
"""

M_earth = 5.9722e24 # [kg]
# Período de um dia em segundos
T_dia_sinodico = 24 * 60 * 60  # [s]
T_dia_sideral = T_dia_sinodico / (1 + 1/365.25)

# Raio da Terra em km
R_terra = 6378.137

omega_GEO = 2 * np.pi / T_dia_sideral
mu = M_earth * Orbit.G
R_GEO = (mu / omega_GEO**2)**(1/3)

# Criar órbita com período de um dia
orbita_GEO = Orbit(m1=M_earth, a = R_GEO, e = 0, body1radius=R_terra)  # Aproximadamente 35786 km de altitude

# Calcular altitude e velocidade
zGEO = orbita_GEO.a - R_terra  # Altitude em km
vGEO = orbita_GEO.v_at_r(orbita_GEO.a)  # Velocidade em km/s

print(f"Altitude da órbita geoestacionária: {zGEO:.3f} km")
print(f"Velocidade na órbita geoestacionária: {vGEO:.3f} km/s")
print(f"Período orbital: {orbita_GEO.T/3600:.3f} horas")

# Plotar a órbita
orbita_GEO.plot()