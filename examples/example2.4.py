from astronautica import Orbit, Body
import numpy as np
import matplotlib.pyplot as plt

"""
Plot the speed v and period T of a satellite in a circular LEO as a function of altitude z.
"""

earth = Body(name="earth")

# Criar array de altitudes de 0 a 2000 km
z = np.linspace(0, 2000, 1000)

# Arrays para armazenar velocidades e períodos
v = np.zeros_like(z)
T = np.zeros_like(z)

# Calcular velocidade e período para cada altitude
for i, altitude in enumerate(z):
    # Criar órbita circular na altitude especificada
    raio = earth.radius + altitude
    orbita = Orbit.from_apsis(earth, rp=raio, ra=raio)
    
    # Armazenar velocidade e período
    v[i] = orbita.v_at_r(raio)
    T[i] = orbita.T / 60  # Convertendo período para minutos

# Criar figura com dois subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Plotar velocidade
ax1.plot(z, v)
ax1.set_xlabel('Altitude (km)')
ax1.set_ylabel('Velocidade (km/s)')
ax1.set_title('Velocidade Orbital vs Altitude')
ax1.grid(True)

# Plotar período
ax2.plot(z, T)
ax2.set_xlabel('Altitude (km)')
ax2.set_ylabel('Período (minutos)')
ax2.set_title('Período Orbital vs Altitude')
ax2.grid(True)

plt.tight_layout()
plt.show()
