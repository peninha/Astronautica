from astronautica import Orbit, Body, Plotter
import numpy as np

"""
Calculate the altitude zGEO and speed vGEO of a geostationary earth satellite.
"""

earth = Body(name="earth")

# Período de um dia em segundos
T_dia_sinodico = 24 * 60 * 60  # [s]
T_dia_sideral = T_dia_sinodico / (1 + 1/365.25)

omega_GEO = 2 * np.pi / T_dia_sideral
mu = earth.mu()
R_GEO = (mu / omega_GEO**2)**(1/3)

# Criar órbita com período de um dia
orbita_GEO = Orbit.from_apsis(earth, rp=R_GEO, ra=R_GEO)  # Aproximadamente 35786 km de altitude

# Calcular altitude e velocidade
zGEO = orbita_GEO.a - earth.radius  # Altitude em km
vGEO = orbita_GEO.v_at_r(orbita_GEO.a)  # Velocidade em km/s

print(f"Altitude da órbita geoestacionária: {zGEO:.3f} km")
print(f"Velocidade na órbita geoestacionária: {vGEO:.3f} km/s")
print(f"Período orbital: {orbita_GEO.T/3600:.3f} horas")

# Plotar a órbita
Plotter(plot3d=True).plot_orbit(orbita_GEO)