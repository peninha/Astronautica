import numpy as np
import matplotlib.pyplot as plt

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
def f(t, y, m1, m2):
    # Problema dos 2 corpos
    G = 6.67430e-20  # Constante gravitacional [km^3*kg^-1*s^-2)]
    R1 = y[0:3]      # Posição da partícula 1
    R2 = y[3:6]      # Posição da partícula 2
    V1 = y[6:9]      # Velocidade da partícula 1
    V2 = y[9:12]     # Velocidade da partícula 2

    # Calcula a distância entre as partículas
    r = np.linalg.norm(R2 - R1)
    r3 = r ** 3

    # Calcula as acelerações devido à gravidade
    A1 = G * m2 * (R2 - R1) / r3  # [km/s^2]
    A2 = G * m1 * (R1 - R2) / r3  # [km/s^2]

    # Retorna a derivada do vetor de estado
    return np.concatenate((V1, V2, A1, A2))

# Método de Runge-Kutta de quarta ordem
def runge_kutta(f, t0, tf, m1, m2, y0, h=0.01):
    t = t0
    y = y0.astype(float)  # Garantir que y seja float para evitar erros de tipo
    resultados = [(t, y.copy())]
    
    n = int((tf - t0) / h)  # Número de passos

    for _ in range(n):
        k1 = h * f(t, y, m1, m2)
        k2 = h * f(t + h / 2, y + k1 / 2, m1, m2)
        k3 = h * f(t + h / 2, y + k2 / 2, m1, m2)
        k4 = h * f(t + h, y + k3, m1, m2)
        
        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        t += h
        resultados.append((t, y.copy()))
    
    return resultados

# Parâmetros iniciais
t0 = 0 # Tempo inicial [s]
tf = 400 # Tempo final [s]

# Partícula 1
m1 = 1e26 # [kg]
R1 = np.array([0, 0, 0]) # [km]
V1 = np.array([10, 20, 30]) # [km/s]

# Partícula 2
m2 = m1 # [kg]
R2 = np.array([3000, 0, 0]) # [km]
V2 = np.array([0, 40, 0]) # [km/s]

# Vetor de estado inicial
y0 = np.concatenate((R1, R2, V1, V2))

h = 0.01 #integration step

# Executando o método de Runge-Kutta
resultados_rk = runge_kutta(f, t0, tf, m1, m2, y0, h)

# Separando os resultados para plotar (Runge-Kutta manual)
t_vals = np.array([t for t, _ in resultados_rk])
R1_vals = np.array([y[0:3] for _, y in resultados_rk])  # Posição da partícula 1 ao longo do tempo
R2_vals = np.array([y[3:6] for _, y in resultados_rk])  # Posição da partícula 2 ao longo do tempo


#%% Plot Referencial Externo
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(R1_vals[:, 0], R1_vals[:, 1], R1_vals[:, 2], label="Partícula 1")
ax.plot(R2_vals[:, 0], R2_vals[:, 1], R2_vals[:, 2], label="Partícula 2")

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


#%% Plot Referencial Centro de Massa
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

CM = (m1*R1_vals + m2*R2_vals)/(m1+m2)

R1_relative_to_CM = R1_vals - CM
R2_relative_to_CM = R2_vals - CM
ax.plot(R1_relative_to_CM[:, 0], R1_relative_to_CM[:, 1], R1_relative_to_CM[:, 2], label="Partícula 1 - Ref CM")
ax.plot(R2_relative_to_CM[:, 0], R2_relative_to_CM[:, 1], R2_relative_to_CM[:, 2], label="Partícula 2 - Ref CM")

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


#%% Plot Referencial R1
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

R2_relative_to_R1 = R2_vals - R1_vals
ax.plot(R2_relative_to_R1[:, 0], R2_relative_to_R1[:, 1], R2_relative_to_R1[:, 2], label="Partícula 2 - Referencial de R1")

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


#%% Plot Referencial R2
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

R1_relative_to_R2 = R1_vals - R2_vals
ax.plot(R2_relative_to_R1[:, 0], R2_relative_to_R1[:, 1], R2_relative_to_R1[:, 2], label="Partícula 1 - Referencial de R2")

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