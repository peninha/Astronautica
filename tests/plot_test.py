from astronautica import Orbit, Body, Plotter, Frame
import numpy as np

earth = Body(name="earth")

e = 0.8
h = 100000

orbita = Orbit.from_elements(earth, e=e, h=h, Omega0=140, i=-30, omega0=90, theta0=-30, t0_clock=100)

# Plotar a órbita
plotter = Plotter(frame="bodycentric", plot3d=True)
plotter.plot_orbit(orbita)

print(orbita)