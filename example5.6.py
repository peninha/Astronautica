from orbitalmechanics import Orbit
import numpy as np

"""
find the local sidereal time (in degrees) of Tokyo, Japan, on March 3, 2004, at 4:30:00 UT. The east
longitude of Tokyo is 139.80°. (This places Tokyo nine time zones ahead of Greenwich, so the local
time is 1:30 in the afternoon.)
"""

theta_g = Orbit.greenwich_sideral_from_date_UT(2004, 3, 3, 4, 30, 0)
theta_l = Orbit.local_sideral_from_theta_g_and_longitude(theta_g, 139.80)

print(f"Local sidereal time: {theta_l} degrees")

