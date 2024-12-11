import matplotlib
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot()
#ax = fig.add_subplot(projection='3d')
ax.plot([1, 2, 3], [4, 5, 6], label='teste')
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.axis('equal')
plt.title("Gráfico em janela externa")
plt.tight_layout()
plt.get_current_fig_manager().window.showMaximized()
plt.show()