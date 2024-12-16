from orbitalmechanics import Orbit
import numpy as np

"""
At the instant the Greenwich sidereal time is θG = 126.7°, the geocentric equatorial
position vector of the International Space Station is:
 r = (-5368, -1784, 3691) km
Find its topocentric right ascension and declination at sea level (H = 0),
latitude ϕ = 20°, and east longitude Λ = 60°.
"""

r = np.array([-5368, -1784, 3691])
theta_g = 126.7
phi = 20
longitude = 60

