from orbitalmechanics import Orbit
import numpy as np

"""
Find the elapsed time between October 4, 1957 UT 19:26:24, and the date of the previous example.
"""

jd1 = Orbit.jd_from_calendar_date(1957, 10, 4, 19, 26, 24)

jd2 = Orbit.jd_from_calendar_date(2004, 5, 12, 14, 45, 30)

print(f"Elapsed time: {(jd2 - jd1) / 365.25} years")

