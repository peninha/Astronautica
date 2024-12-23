from astronautica import Orbit, Body, Plotter
import numpy as np
from scipy.optimize import fsolve

"""
Let a satellite be in a 500 km by 5000 km orbit with its apse line parallel to the line from the earth to the sun,
as shown in Fig. 3.9. Find the time that the satellite is in the earth’s shadow if
(a) the apogee is toward the sun
(b) the perigee is toward the sun.
"""

earth = Body(name="earth")
rp = 500 + earth.radius # [km]
ra = 5000 + earth.radius # [km]

orbita = Orbit.from_elements(earth, rp=rp, ra=ra, theta0=0)

"""
Condição para sombra:
r * sen(theta) = R_earth
Da equação da órbita: r = h**2/mu * 1/(1 + e*cos(theta))
Temos: h**2/mu * sen(theta) - R_earth * (1 + e * cos(theta)) = 0
"""
h = orbita.h
mu = orbita.mu
e = orbita.e
R_earth = earth.radius

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
t1 = orbita.t_clock_at_theta(theta1)*2
print(f"Tempo de sombra: {t1/60:.4f} minutos")

print("b) perigeu para o sol")
t2 = orbita.T - orbita.t_clock_at_theta(theta2)*2
print(f"Tempo de sombra: {t2/60:.4f} minutos")

orbita.add_orbital_position(theta=theta1, name="sombra com apogeu para o sol") 
orbita.add_orbital_position(theta=theta2, name="sombra com perigeu para o sol")

plot = Plotter(frame="bodycentric", plot3d=True)
plot.plot_orbit(orbit=orbita)