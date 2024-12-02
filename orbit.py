import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def set_axes_equal(ax: plt.Axes):
    """Set 3D plot axes to equal scale.

    Make axes of 3D plot have equal scale so that spheres appear as
    spheres and cubes as cubes.  Required since `ax.axis('equal')`
    and `ax.set_aspect('equal')` don't work on 3D.
    """
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)

def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])

# Defina a função f(t, y), que é a derivada dy/dt
def f(t, y, mu):
    # orbita
    r = y[0:3]      # Posição do corpo em órbita
    v = y[3:6]      # Velocidade do corpo em órbita

    # Calcula a distância entre as partículas
    r_norm = np.linalg.norm(r)
    r_norm3 = r_norm**3

    # Calcula as acelerações devido à gravidade
    a = -mu * r / r_norm3  # [km/s^2]

    # Retorna a derivada do vetor de estado
    return np.concatenate((v, a))

# Parâmetros iniciais
G = 6.67430e-20  # Constante gravitacional [km^3*kg^-1*s^-2)]
t0 = 0 # Tempo inicial [s]
tf = 14400 # Tempo final [s]

# Corpo de referência
M = 5.9722e24 # Corpo de Referência (Earth) [kg]

# Corpo em órbita
m = 1000 # [kg]
r = np.array([8000, 0, 6000]) # [km]
v = np.array([0, 7, 0]) # [km/s]

# Vetor de estado inicial
y0 = np.concatenate((r, v))

h = 1 #time resolution [s]
mu = G*(M+m)

# Executando o método de Runge-Kutta
#resultados_rk = runge_kutta(f, t0, tf, m1, m2, y0, h)

# Executando o método de integração usando solve_ivp
sol = solve_ivp(f, [t0, tf], y0, args=(mu,), method='RK45', t_eval=np.linspace(t0, tf, int((tf - t0) / h)), rtol=1e-9, atol=1e-12)
#sol = solve_ivp(f, [t0, tf], y0, args=(m1, m2), method='DOP853', t_eval=np.linspace(t0, tf, int((tf - t0) / h)), rtol=1e-9, atol=1e-12)
#sol = solve_ivp(f, [t0, tf], y0, args=(m1, m2), method='Radau', t_eval=np.linspace(t0, tf, int((tf - t0) / h)), rtol=1e-9, atol=1e-12)

# Separando os resultados para plotar (solve_ivp)
R1_vals = sol.y[0:3, :].T  # Posição do córpo em órbita ao longo do tempo

#%% Plot Referencial M
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plotando a trajetória do corpo em órbita
ax.plot(R1_vals[:, 0], R1_vals[:, 1], R1_vals[:, 2], label="Corpo em Órbita")

# Plotando uma esfera simplificada representando a Terra
# Usando menos pontos para uma representacao mais leve
u = np.linspace(0, 2 * np.pi, 30)
v = np.linspace(0, np.pi, 15)
x = 6371 * np.outer(np.cos(u), np.sin(v))  # Raio da Terra [km]
y = 6371 * np.outer(np.sin(u), np.sin(v))
z = 6371 * np.outer(np.ones(np.size(u)), np.cos(v))

ax.plot_wireframe(x, y, z, color='b', alpha=0.5)

# Plotando os eixos X, Y, Z
ax.quiver(0, 0, 0, 10000, 0, 0, color='r', arrow_length_ratio=0.05, label='Eixo X')
ax.quiver(0, 0, 0, 0, 10000, 0, color='g', arrow_length_ratio=0.05, label='Eixo Y')
ax.quiver(0, 0, 0, 0, 0, 10000, color='k', arrow_length_ratio=0.05, label='Eixo Z')

ax.legend()

# Manter a mesma escala em todos os eixos
ax.set_box_aspect([1,1,1]) # IMPORTANT - this is the new, key line
set_axes_equal(ax) # IMPORTANT - this is also required

# Removendo os eixos para uma aparência mais clean
ax.set_axis_off()

plt.show()