from orbitalmechanics import Orbit
import numpy as np

"""
If the direction cosine matrix for the transformation from xyz to x'y'z' is
Q = 0.64050 0.75309 -0.15038
    0.76737 -0.63530 0.086823
    -0.30152 -0.17101 -0.98481
find the angles α, β, and γ of the classical Euler sequence.
"""

Q = np.array([[0.64050, 0.75309, -0.15038],
              [0.76737, -0.63530, 0.086823],
              [-0.030152, -0.17101, -0.98481]])

alpha, beta, gamma = Orbit.Euler_angles_from_Q(Q, pattern="yaw-pitch-roll")
print(alpha, beta, gamma)
