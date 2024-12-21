from orbitalmechanics import Orbit
import numpy as np
from scipy.optimize import fsolve

"""
Let a satellite be in a 500 km by 5000 km orbit with its apse line parallel to the line from the earth to the sun,
as shown in Fig. 3.9. Find the time that the satellite is in the earth’s shadow if
(a) the apogee is toward the sun
(b) the perigee is toward the sun.
"""

M_earth = 5.974e24 # [kg]
Earth_radius = 6378 # [km]
rp = 500 + Earth_radius # [km]
ra = 5000 + Earth_radius # [km]

orbita = Orbit(m1=M_earth, rp=rp, ra=ra, body1radius=Earth_radius)

"""
Condição para sombra:
r * sen(theta) = R_earth
Da equação da órbita: r = h**2/mu * 1/(1 + e*cos(theta))
Temos: h**2/mu * sen(theta) - R_earth * (1 + e * cos(theta)) = 0
"""
h = orbita.h
mu = orbita.mu
e = orbita.e
R_earth = Earth_radius

def f(theta):
    return h**2/mu * np.sin(theta) - R_earth*(1 + e*np.cos(theta))

# Procurar a primeira raiz (entrada na sombra)
theta1 = fsolve(f, np.pi/2)[0]  # chute inicial próximo a π/2
theta1 = np.mod(np.degrees(theta1), 360)

# Procurar a segunda raiz (saída da sombra)
theta2 = fsolve(f, 3*np.pi/2)[0]  # chute inicial próximo a 3π/2
theta2 = np.mod(np.degrees(theta2), 360)

print(f"Ângulos de entrada/saída da sombra: {theta1:.4f}, {theta2:.4f} graus")

print("a) apogeu para o sol")
t = orbita.t_at_theta(theta1)*2
print(f"Tempo de sombra: {t/60:.4f} minutos")

print("b) perigeu para o sol")
t = orbita.T - orbita.t_at_theta(theta2)*2
print(f"Tempo de sombra: {t/60:.4f} minutos")

orbita.add_orbital_position(theta1, name="sombra com apogeu para o sol") 
orbita.add_orbital_position(theta2, name="sombra com perigeu para o sol")
orbita.plot(plot_positions=True)
