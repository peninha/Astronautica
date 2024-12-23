from astronautica import Orbit, Body, Plotter
import numpy as np

"""
An earth orbit has an eccentricity of 0.3, an angular momentum of 60000 km2/s, 
and a true anomaly of 120°. What are the position vector r and velocity vector v
in the perifocal frame of reference?
"""
e = 0.3
h = 60000
theta = 120

earth = Body(name="earth")  

orbit = Orbit.from_elements(earth, e=0.3, h=60000, theta0=theta)

state_vectors = orbit.state_vectors_at_theta(theta, "perifocal")
print(state_vectors)

Plotter(frame="perifocal", plot3d=False).plot_orbit(orbit=orbit)