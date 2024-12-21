from orbitalmechanics import Orbit
import numpy as np

"""
The geocentric position vectors of a space object at three successive times are
r1 = (-294.32, 4265.1, 5986.7) km
r2 = (-1365.5, 3637.6, 6346.8) km
r3 = (-2940.3, 2473.7, 6555.8) km
Determine the classical orbital elements using Gibbs method.
"""
M_earth = 5.9722e24 # [kg]
R_terra = 6378.137

r1 = np.array([-294.32, 4265.1, 5986.7])
r2 = np.array([-1365.5, 3637.6, 6346.8])
r3 = np.array([-2940.3, 2473.7, 6555.8])

orbita = Orbit.init_from_3_vectors(r1, r2, r3, m1=M_earth, m2=0, body1radius=R_terra)

print(orbita)

orbita.plot(frame="bodycentric", points=True, velocities=True, positions=True, trajectory=False, plot3d=True)

