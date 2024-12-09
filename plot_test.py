import matplotlib
#matplotlib.use('qt5agg')
print(matplotlib.get_configdir())
print(matplotlib.get_backend())
import matplotlib.pyplot as plt

print("A")

plt.plot([1, 2, 3], [4, 5, 6])
plt.title("Gráfico em janela externa")
plt.show()

print("B")
