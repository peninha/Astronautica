from astronautica import Orbit, Body, Plotter
import numpy as np

"""
The geocentric position vectors of a space object at three successive times are
r1 = (-294.32, 4265.1, 5986.7) km
r2 = (-1365.5, 3637.6, 6346.8) km
r3 = (-2940.3, 2473.7, 6555.8) km
Determine the classical orbital elements using Gibbs method.
"""
earth = Body("earth")

r1 = np.array([-294.32, 4265.1, 5986.7])
r2 = np.array([-1365.5, 3637.6, 6346.8])
r3 = np.array([-2940.3, 2473.7, 6555.8])

orbita = Orbit.from_3_vectors(earth, r1, r2, r3)

#print(orbita)

plotter = Plotter(plot3d=True)
plotter.plot_orbit(orbita, frame="bodycentric", points=True, velocities=True, positions=True)

