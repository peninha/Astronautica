import numpy as np
import matplotlib.pyplot as plt

# Defina a função f(t, y), que é a derivada dy/dt
def f(t, y):
    return -2 * t * y  # Exemplo de função

# Método de Runge-Kutta de quarta ordem
def runge_kutta(f, t0, y0, h, n):
    t = t0
    y = y0
    resultados = [(t, y)]
    
    for _ in range(n):
        k1 = h * f(t, y)
        k2 = h * f(t + h/2, y + k1/2)
        k3 = h * f(t + h/2, y + k2/2)
        k4 = h * f(t + h, y + k3)
        
        y += (k1 + 2*k2 + 2*k3 + k4) / 6
        t += h
        resultados.append((t, y))
    
    return np.array(resultados)

# Parâmetros
t0 = 0          # Tempo inicial
y0 = 1          # Valor inicial de y
h = 0.1         # Passo de integração
n = 100         # Número de passos

# Executando o método
resultados = runge_kutta(f, t0, y0, h, n)

# Plotando os resultados
t_vals = resultados[:, 0]
y_vals = resultados[:, 1]

plt.plot(t_vals, y_vals, label="RK4")
plt.xlabel("t")
plt.ylabel("y")
plt.title("Método de Runge-Kutta 4ª Ordem")
plt.legend()
plt.grid(True)
plt.show()
