from astronautica import Orbit, Body, Plotter
import numpy as np

"""
An earth satellite has the following position and velocity vectors at a given instant:
r = 7000p_vec + 9000q_vec km
v = -3.3472p_vec + 9.1251q_vec km/s
Calculate the specific angular momentum h, the true anomaly θ, and the eccentricity e.
"""
earth = Body(name="earth")

# Vetores de posição e velocidade no sistema perifocal
r_vec = [7000, 9000, 0]  # km
v_vec = [-3.3472, 9.1251, 0]  # km/s

# Criar órbita a partir dos vetores
orbita = Orbit.from_state_vectors(earth, r_vec_bc=r_vec, v_vec_bc=v_vec)

# Calcular momento angular específico
print(f"\nMomento angular específico: {orbita.h:.2f} km²/s")

# Calcular anomalia verdadeira
r = np.linalg.norm(r_vec)  # magnitude do vetor posição
theta = orbita.theta_at_r(r)
print(f"Anomalia verdadeira: {theta[0]:.2f}°")

# Calcular excentricidade
print(f"Excentricidade: {orbita.e:.4f}")

# Plotar a órbita com o ponto atual
Plotter(frame="bodycentric", plot3d=True).plot_orbit(orbit=orbita)
