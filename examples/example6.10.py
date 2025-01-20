from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
Determine the required launch azimuth for the sun-synchronous satellite of Example 4.9 if it is
launched from Vandenburgh AFB on the California coast (latitude = 34.5°N).
"""
earth = Body("spherical_earth")

latitude = 34.5
inclination = 98.43

launch_azimuth_prograde = Orbit.launch_azimuth(latitude, inclination, to_north=True)
launch_azimuth_retrograde = Orbit.launch_azimuth(latitude, inclination, to_north=False)

print(f"Launch azimuth (to north): {launch_azimuth_prograde:.2f}°")
print(f"Launch azimuth (to south): {launch_azimuth_retrograde:.2f}°")

