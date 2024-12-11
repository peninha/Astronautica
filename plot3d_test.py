import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(projection='3d')

ax.plot([1, 2, 3], [4, 5, 6], [7, 8, 9], label='teste')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 1 + 2 * np.outer(np.cos(u), np.sin(v))
y = 4 + 2 * np.outer(np.sin(u), np.sin(v))
z = 7 + 2 * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='red', alpha=0.4, label='bola')

ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')

ax.axis('equal')
#ax.grid(False)
plt.title("Gráfico em janela externa")
plt.tight_layout()
plt.get_current_fig_manager().window.showMaximized()
plt.show()