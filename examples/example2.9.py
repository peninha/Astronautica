from astronautica import Orbit, Body, Plotter
import numpy as np

"""
The perigee radius of a satellite in a parabolic geocentric trajectory of Fig. 2.24 is 7000 km. Find the distance d between
points P1 and P2 on the orbit, which are 8000 km and 16000 km, respectively, from the center of the earth.
"""

earth = Body(name="earth")

# Criar órbita parabólica (e = 1) com raio do perigeu dado
rp = 7000  # km
orbita = Orbit.from_elements(earth, rp=rp, e=1)

# Pontos dados
r1 = 8000  # km
r2 = 16000  # km

# Calcular anomalias verdadeiras para cada ponto
theta1 = orbita.theta_at_r(r1)[0]  # Pegamos apenas o primeiro valor pois é simétrico
theta2 = orbita.theta_at_r(r2)[0]

# Converter para coordenadas cartesianas
x1 = r1 * np.cos(np.radians(theta1))
y1 = r1 * np.sin(np.radians(theta1))

x2 = r2 * np.cos(np.radians(theta2))
y2 = r2 * np.sin(np.radians(theta2))

# Calcular distância entre os pontos
d = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

print(orbita)

print(f"Distância entre P1 e P2: {d:.2f} km")

# Plotar a órbita com os pontos
orbita.add_orbital_position(theta=theta1, name='P1')
orbita.add_orbital_position(theta=theta2, name='P2')
Plotter(frame="bodycentric", plot3d=True).plot_orbit(orbita)
