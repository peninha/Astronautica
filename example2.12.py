from orbitalmechanics import Orbit
import numpy as np

"""
An earth satellite has the following position and velocity vectors at a given instant:
r = 7000p_vec + 9000q_vec km
v = -3.3472p_vec + 9.1251q_vec km/s
Calculate the specific angular momentum h, the true anomaly θ, and the eccentricity e.
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137 # [km]
# Vetores de posição e velocidade no sistema perifocal
r_vec = [7000, 9000, 0]  # km
v_vec = [-3.3472, 9.1251, 0]  # km/s

# Criar órbita a partir dos vetores
orbita = Orbit.init_from_state_vectors(m1=M_earth, r_vec=r_vec, v_vec=v_vec, body1radius=R_terra)

# Calcular momento angular específico
print(f"\nMomento angular específico: {orbita.h:.2f} km²/s")

# Calcular anomalia verdadeira
r = np.linalg.norm(r_vec)  # magnitude do vetor posição
theta = orbita.theta_at_r(r)
print(f"Anomalia verdadeira: {theta[0]:.2f}°")

# Calcular excentricidade
print(f"Excentricidade: {orbita.e:.4f}")

# Plotar a órbita com o ponto atual
orbita.plot(plot_velocities=True)
