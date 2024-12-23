import numpy as np

"""
In Fig. 4.10, the x' axis is defined by the line segment O'P.
The x'y' plane is defined by the intersecting line segments O'P and O'Q.
The z' axis is normal to the plane of O'P and O'Q and obtained by rotating O'P toward O'Q and using the righthand rule.
O' = [3, 1, 2]
P = [-5, 5, 4]
Q = [-6, 3, 5]

(a) Find the direction cosine matrix [Q].
(b) If {v} = [2, 4, 6]T, find {v'}.
(c) If {v'} = [2, 4, 6]T, find {v}.
"""
O_prime = np.array([3, 1, 2])
P = np.array([-5, 5, 4])
Q = np.array([-6, 3, 5])

O_prime_P = P - O_prime
O_prime_Q = Q - O_prime

i_prime = O_prime_P/np.linalg.norm(O_prime_P)
k_prime = np.cross(O_prime_P, O_prime_Q)/np.linalg.norm(np.cross(O_prime_P, O_prime_Q))
j_prime = np.cross(k_prime, i_prime)

Q_prime = np.array([i_prime, j_prime, k_prime])

print(Q_prime)

v = np.array([2, 4, 6])
v_prime = Q_prime @ v
print(v_prime)

v_prime = np.array([2, 4, 6])
v = Q_prime.T @ v_prime
print(v)

