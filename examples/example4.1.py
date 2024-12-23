from astronautica import Frame
import numpy as np

"""
If the position vector of the International Space Station in the geocentric equatorial frame is
r = -5368I - 1784J + 3691K (km) , what are its right ascension and declination?
"""

r_vec = np.array([-5368, -1784, 3691]) # [km]
ra, dec = Frame.convert_cartesian_to_ra_dec(r_vec)
print(f"Right ascension: {ra:.4f}°")
print(f"Declination: {dec:.4f}°")
