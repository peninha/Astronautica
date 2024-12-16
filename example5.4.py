from orbitalmechanics import Orbit
import numpy as np

"""
What is the Julian day number for May 12, 2004, at 14:45:30 UT?
"""

jd = Orbit.jd_from_calendar_date(2004, 5, 12, 14, 45, 30)

print(jd)

