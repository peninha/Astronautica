import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar, fsolve

class OrbitalBase:
    """
    Base class for all orbital mechanics classes.
    """
    G = 6.67430e-20  # Gravitational constant in km^3/kg/s^2

    def __init__(self):
        self.points_in_orbital_plane = []
        self.positions = []
        self.trajectory_points = []
    
    ########################   Calculations   ########################
    #####   mu   #####
    @classmethod
    def mu_from_m1_m2(cls, m1, m2=0):
        """
        Calculates the gravitational parameter mu based on the masses of two bodies.
        
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body, default is 0 (kg)
        """
        return cls.G * (m1 + m2)
    
    #####   r, theta   #####
    @classmethod
    def r_theta_from_state_vectors(cls, r_vec, v_vec, mu=None, m1=None, m2=0):
        """
        Calculates the distance and true anomaly from the position and velocity vectors,
        considering any reference frame centered at the primary body.

        :param r_vec: Position vector (km)
        :param v_vec: Velocity vector (km/s)
        :param mu: Gravitational parameter (km^3/s^2), default is None
        :param m1: Mass of the primary body (kg), default is None
        :param m2: Mass of the secondary body (kg), default is 0
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Provide mu or m1 (and optionally m2) to calculate mu.")

        orb_params = cls.orbital_parameters_from_state_vectors(r_vec, v_vec, mu, m1, m2)
        return orb_params["r"], orb_params["theta"]
    
    #####   Orbital parameters   #####
    @classmethod
    def orbital_parameters_from_state_vectors(cls, r_vec, v_vec, mu=None, m1=None, m2=0):
        """
        Calculates the orbital parameters from the position and velocity vectors,
        considering any reference frame centered at the primary body.
        
        :param r_vec: Position vector (km)
        :param v_vec: Velocity vector (km/s)
        :param mu: Gravitational parameter (km^3/s^2), default is None
        :param m1: Mass of the primary body (kg), default is None
        :param m2: Mass of the secondary body (kg), default is 0
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Provide mu or m1 (and optionally m2) to calculate mu.")

        # Angular momentum vector
        h_vec = cls.h_vec_from_mu_state_vectors(mu, r_vec, v_vec)
        h = np.linalg.norm(h_vec)

        # Inclination   
        i = cls.i_from_mu_h_vec(mu, h_vec)

        # Node vector
        n_vec = cls.n_vec_from_mu_h_vec(mu, h_vec)
        n = np.linalg.norm(n_vec)

        # Right ascension of the ascending node
        Omega = cls.Omega_from_mu_n_vec(mu, n_vec)

        # Eccentricity vector
        e_vec = cls.e_vec_from_mu_state_vectors(mu, r_vec, v_vec)
        e = np.linalg.norm(e_vec)

        # Argument of periapsis
        omega = cls.omega_from_mu_e_vec_n_vec(mu, e_vec, n_vec)

        # True anomaly
        theta = cls.theta_from_mu_h_vec_e_vec_r_vec(mu, h_vec, e_vec, r_vec)

        # Distance
        r = np.linalg.norm(r_vec)

        # Velocity
        v = np.linalg.norm(v_vec)

        # Rotation matrix
        Q_bc_peri = cls.Q_from_Euler_angles(Omega, i, omega, pattern="classical")

        # Other parameters
        rp = cls.rp_from_mu_h_e(mu, h, e)
        ra = cls.ra_from_mu_h_e(mu, h, e)
        a = cls.a_from_rp_ra(rp, ra)
        epsilon = cls.epsilon_from_mu_a(mu, a)
        p = cls.p_from_mu_h(mu, h)
        b = cls.b_from_a_e(a, e)
        T = cls.T_from_mu_a(mu, a)
        alpha = cls.alpha_from_a(a)

        orb_params = {
            "mu": mu,
            "m1": m1,
            "m2": m2,
            "h": h,
            "e": e,
            "i": i,
            "Omega": Omega,
            "omega": omega,
            "theta": theta,
            "r_vec": r_vec,
            "v_vec": v_vec,
            "h_vec": h_vec,
            "e_vec": e_vec,
            "n_vec": n_vec,
            "r": r,
            "v": v,
            "n": n,
            "Q_bc_peri": Q_bc_peri,
            "rp": rp,
            "ra": ra,
            "a": a,
            "b": b,
            "p": p,
            "T": T,
            "epsilon": epsilon,
            "alpha": alpha,
        }
        return orb_params

    #####   i   #####
    @staticmethod
    def i_from_mu_h_vec(mu, h_vec):
        """
        Calculates the inclination using gravitational parameter and specific angular momentum vector.
        """
        return np.degrees(np.arccos(h_vec[2]/np.linalg.norm(h_vec)))

    #####   N   #####
    @staticmethod
    def n_vec_from_mu_h_vec(mu, h_vec):
        """
        Calculates the node vector using gravitational parameter and specific angular momentum vector.
        """
        return np.cross(np.array([0, 0, 1]), h_vec)

    #####   Omega   #####
    @staticmethod
    def Omega_from_mu_n_vec(mu, n_vec):
        """
        Calculates the right ascension of the ascending node using gravitational parameter and node vector.
        """
        n = np.linalg.norm(n_vec)
        if n == 0:
            Omega = 0
        else:
            cos_Omega = n_vec[0]/n
            sin_Omega = n_vec[1]/n
            Omega = np.mod(np.degrees(np.arctan2(sin_Omega, cos_Omega)), 360)
        return Omega

    #####   omega   #####
    @staticmethod
    def omega_from_mu_e_vec_n_vec(mu, e_vec, n_vec):
        """
        Calculates the argument of periapsis using gravitational parameter, eccentricity vector and node vector.
        """
        e = np.linalg.norm(e_vec)
        n = np.linalg.norm(n_vec)
        if n == 0 or e == 0:
            omega = 0
        else:
            omega = np.degrees(np.arccos(np.dot(n_vec, e_vec)/(n * e))) 
            if e_vec[2] < 0:
                omega = 360 - omega
        return omega

    #####   theta   #####
    @staticmethod
    def theta_from_mu_h_vec_e_vec_r_vec(mu, h_vec, e_vec, r_vec):
        """
        Calculates the true anomaly using gravitational parameter, eccentricity vector and position vector.
        """
        h = np.linalg.norm(h_vec)
        e = np.linalg.norm(e_vec)
        r = np.linalg.norm(r_vec)
        if e == 0:
            theta = 0
        else:
            cos_theta = np.dot(e_vec, r_vec)/(e * r)
            sin_theta = np.dot(h_vec, np.cross(e_vec, r_vec))/(h * e * r)
            theta = np.mod(np.degrees(np.arctan2(sin_theta, cos_theta)), 360)
        return theta

    #####   rp   #####
    @staticmethod
    def rp_from_a_e(a, e):
        """
        Calculates the periapsis distance using semi-major axis and eccentricity.
        Note: For hyperbolic orbits (e > 1), 'a' must be a negative value.
        """
        return a * (1 - e)

    @staticmethod
    def rp_from_mu_h_e(mu, h, e):
        """
        Calculates the periapsis distance using gravitational parameter, specific angular momentum and eccentricity.
        """
        return h**2 / (mu * (1 + e))

    #####   ra   #####
    @staticmethod
    def ra_from_a_e(a, e):
        """
        Calculates the apoapsis distance using semi-major axis and eccentricity.
        Note: For hyperbolic orbits (e > 1), 'a' must be a negative value, 
        which will result in a negative distance.
        """
        if e == 1:
            return float('inf')
        return a * (1 + e)
    
    @staticmethod
    def ra_from_mu_h_e(mu, h, e):
        """
        Calculates the apoapsis distance using gravitational parameter, specific angular momentum and eccentricity.
        """
        if e == 1:
            return float('inf')
        else:
            return h**2 / (mu * (1 - e))

    #####   a   #####
    @staticmethod
    def a_from_rp_ra(rp, ra):
        """
        Calculates the semi-major axis using periapsis and apoapsis distances.
        Note: For hyperbolic orbits (e > 1), 'ra' must be a negative value, 
        which will result in a negative distance.
        """
        if ra == float('inf'):
            return float('inf')
        return (rp + ra) / 2

    @staticmethod
    def a_from_mu_h_e(mu, h, e):
        """
        Calculates the semi-major axis using gravitational parameter, specific angular momentum and eccentricity.
        """
        if e == 1:
            return float('inf')
        return h**2 / (mu * (1 - e**2))
    
    @staticmethod
    def a_from_mu_r_v(mu, r, v):
        """
        Calculates the semi-major axis using gravitational parameter, distance and velocity.
        """
        return 1 / (2/r - v**2/mu)

    #####   e   #####
    @staticmethod
    def e_from_rp_ra(rp, ra):
        """
        Calculates the eccentricity using periapsis and apoapsis distances.
        Note: For hyperbolic orbits (e > 1), 'ra' must be a negative value. 
        """
        if ra == float('inf'):
            return 1
        return (ra - rp) / (ra + rp)
    
    @staticmethod
    def e_from_mu_h_rp(mu, h, rp):
        """
        Calculates the eccentricity using gravitational parameter, specific angular momentum and periapsis distance.
        """
        return h**2 / (mu * rp) - 1
    
    @staticmethod
    def e_from_mu_r_v_theta(mu, r, v, theta):
        """
        Calculates the eccentricity using gravitational parameter, distance, velocity and true anomaly.
        """
        # solve second degree equation for eccentricity
        u = r*v**2/mu
        c1 = 1
        c2 = np.cos(theta)*(2 - u)
        c3 = 1 - u
        delta = c2**2 - 4*c1*c3
        e = (-c2 + np.sqrt(delta))/(2*c1) # positive root
        return e

    @staticmethod
    def e_vec_from_mu_state_vectors(mu, r_vec, v_vec):
        """
        Calculates the eccentricity vector using gravitational parameter and state vectors.
        """
        r = np.linalg.norm(r_vec)
        v = np.linalg.norm(v_vec)
        vr = np.dot(r_vec, v_vec)/r
        return 1/mu * ( (v**2 - mu/r)*r_vec - r*vr*v_vec )
    

    #####   h   #####
    @staticmethod
    def h_from_mu_a_e(mu, a, e):
        """
        Calculates the specific angular momentum using gravitational parameter, semi-major axis and eccentricity.
        For hyperbolic orbits (e > 1), a should be negative.
        """
        if e == 1:
            raise ValueError("Can't calculate h from 'a' and 'e' for parabolic (e=1) orbits.")
        return np.sqrt(mu * a * (1 - e**2))
    
    @staticmethod
    def h_from_mu_rp_e(mu, rp, e):
        """
        Calculates the specific angular momentum using gravitational parameter, periapsis distance and eccentricity.
        """
        return np.sqrt(mu * rp * (1 + e))

    @staticmethod
    def h_from_mu_r_e_theta(mu, r, e, theta):
        """
        Calculates the specific angular momentum using gravitational parameter, distance, eccentricity and true anomaly.
        """
        return np.sqrt(mu * r * (1 + e*np.cos(theta)))
    
    @staticmethod
    def h_from_mu_state_vectors(mu, r_vec, v_vec):
        """
        Calculates the specific angular momentum using gravitational parameter, position and velocity vectors.
        """
        return np.linalg.norm(OrbitalBase.h_vec_from_mu_state_vectors(mu, r_vec, v_vec))
    
    @staticmethod
    def h_vec_from_mu_state_vectors(mu, r_vec, v_vec):
        """
        Calculates the specific angular momentum vector using gravitational parameter, position and velocity vectors.
        """
        return np.cross(r_vec, v_vec)

    #####   b   #####
    @staticmethod
    def b_from_a_e(a, e):
        """
        Calculates the semi-minor axis using semi-major axis and eccentricity.
        Note: For hyperbolic orbits (e > 1), 'a' must be a negative value.
        """
        if e < 1:
            return a * np.sqrt(1 - e**2)
        elif e == 1:
            return None
        else: # e > 1
            return -a * np.sqrt(e**2 - 1)
    
    #####   p   #####
    @staticmethod
    def p_from_mu_h(mu, h):
        """
        Calculates the semi-latus rectum using gravitational parameter and specific angular momentum.
        """
        return h**2 / mu

    #####   epsilon   #####
    @staticmethod
    def epsilon_from_mu_h_e(mu, h, e):
        """
        Calculates the specific orbital energy using gravitational parameter,
        specific angular momentum and eccentricity.
        """
        return -mu**2 / (2 * h**2) * (1 - e**2)
    
    @staticmethod
    def epsilon_from_mu_a(mu, a):
        """
        Calculates the specific orbital energy using gravitational parameter and semi-major axis.
        Note: For hyperbolic orbits (e > 1), 'a' must be a negative value.
        """
        return -mu / (2 * a)

    #####   T   #####
    @staticmethod
    def T_from_mu_a(mu, a):
        """
        Calculates the orbital period using gravitational parameter and semi-major axis.
        Note: For hyperbolic orbits (e > 1), 'a' must be a negative value.
        """
        if a == float('inf') or a < 0:
            return float('inf')
        return 2 * np.pi * np.sqrt(a**3 / mu)

    #####   alpha   #####
    @staticmethod
    def alpha_from_mu_r_v(mu, r, v):
        """
        Calculates the alpha (inverse of the semi-major axis) parameter using gravitational parameter, distance and velocity.
        """
        return 2/r - v**2/mu

    @staticmethod
    def alpha_from_a(a):
        """
        Calculates the alpha (inverse of the semi-major axis) parameter using semi-major axis.
        """
        if a == float('inf'):
            return 0
        return 1/a

    #####   Stumpff functions   #####
    @staticmethod
    def S(z):
        """
        Stumpff function S(z).
        """
        if z > 0:
            return (np.sqrt(z) - np.sin(np.sqrt(z)))/np.sqrt(z)**3
        elif z < 0:
            return (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/np.sqrt(-z)**3
        else: # z = 0
            return 1/6

    @staticmethod
    def C(z):
        """
        Stumpff function C(z).
        """
        if z > 0:
            return (1 - np.cos(np.sqrt(z)))/z
        elif z < 0:
            return (np.cosh(np.sqrt(-z)) - 1)/(-z)
        else: # z = 0
            return 1/2

    #####   Conversions   #####
    @staticmethod
    def convert_cartesian_to_polar(r_vec_xyz):
        """
        Converts Cartesian coordinates to polar coordinates.
        """
        r = np.sqrt(r_vec_xyz[0]**2 + r_vec_xyz[1]**2)
        theta = np.mod(np.degrees(np.arctan2(r_vec_xyz[1], r_vec_xyz[0])), 360)
        return r, theta

    @staticmethod
    def convert_polar_to_cartesian(r, theta, z=0):
        """
        Converts polar coordinates to Cartesian coordinates.
        """
        theta = np.radians(theta)
        return (r * np.cos(theta), r * np.sin(theta), z)

    @staticmethod
    def convert_cartesian_to_ra_dec(r_vec_xyz):
        """
        Calculates the right ascension and declination from the Cartesian position vector.
        """
        r = np.linalg.norm(r_vec_xyz)
        ra = np.mod(np.degrees(np.arctan2(r_vec_xyz[1], r_vec_xyz[0])), 360)
        dec = np.mod(np.degrees(np.arcsin(r_vec_xyz[2]/r)), 360)
        return ra, dec

    @staticmethod
    def convert_bc_to_rbc(r_vec_bc, omega_body):
        """
        Converts a body-centric vector to a rotating body-centric vector.

        :param r_vec_bc: Body-centric position vector (km)
        :param omega_body: Body rotation rate (rad/s)
        :return: Rotating body-centric position vector (km)
        """
        Q_bc_rbc = OrbitalBase.Q_from_Euler_angles(omega_body, pattern="z")
        return Q_bc_rbc @ r_vec_bc

    #####   Rotation   #####
    @staticmethod
    def Euler_angles_from_Q(Q, pattern="classical"):
        """
        Calculates the Euler angles from the direction cosine matrix.
        """
        if pattern == "313" or pattern == "classical":
            alpha = np.mod(np.degrees(np.arctan2(Q[2, 0], -Q[2, 1])), 360)
            beta = np.degrees(np.arccos(Q[2, 2]))
            gamma = np.mod(np.degrees(np.arctan2(Q[0, 2], Q[1, 2])), 360)
        elif pattern == "321" or pattern == "yaw-pitch-roll":
            alpha = np.mod(np.degrees(np.arctan2(Q[0, 1], Q[0, 0])), 360)
            beta = np.degrees(np.arcsin(-Q[0, 2]))
            gamma = np.mod(np.degrees(np.arctan2(Q[1, 2], Q[2, 2])), 360)
        return alpha, beta, gamma

    @staticmethod
    def Q_from_Euler_angles(alpha, beta=0, gamma=0, pattern="classical"):
        """
        Calculates the direction cosine matrix from the Euler angles.

        """
        alpha = np.radians(alpha)
        beta = np.radians(beta)
        gamma = np.radians(gamma)
        if pattern == "1" or pattern == "x":
            Q = np.array([[1, 0, 0],
                           [0, np.cos(alpha), np.sin(alpha)],
                           [0, -np.sin(alpha), np.cos(alpha)]])
        elif pattern == "2" or pattern == "y":
            Q = np.array([[np.cos(alpha), 0, -np.sin(alpha)],
                           [0, 1, 0],
                           [np.sin(alpha), 0, np.cos(alpha)]])
        elif pattern == "3" or pattern == "z":
            Q = np.array([[np.cos(alpha), np.sin(alpha), 0],
                           [-np.sin(alpha), np.cos(alpha), 0],
                           [0, 0, 1]])
        elif pattern == "313" or pattern == "classical":
            Q = np.array([[np.cos(alpha)*np.cos(gamma) - np.sin(alpha)*np.cos(beta)*np.sin(gamma), np.sin(alpha)*np.cos(gamma) + np.cos(alpha)*np.cos(beta)*np.sin(gamma), np.sin(beta)*np.sin(gamma)],
                          [-np.cos(alpha)*np.sin(gamma) - np.sin(alpha)*np.cos(beta)*np.cos(gamma), -np.sin(alpha)*np.sin(gamma) + np.cos(alpha)*np.cos(beta)*np.cos(gamma), np.sin(beta)*np.cos(gamma)],
                          [np.sin(alpha)*np.sin(beta), -np.cos(alpha)*np.sin(beta), np.cos(beta)]])
        return Q

    ########################   Self methods   ########################

    #######   Miscellaneous   #######
    def type(self):
        """
        Classifies the type of orbit based on eccentricity.
        :return: Type of orbit as a string ("Circular", "Elliptical", "Parabolic", "Hyperbolic")
        """
        if abs(self.e) < 1e-10:  # very close to zero
            return "Circular"
        elif 0 < self.e < 1:
            return "Elliptical"
        elif self.e == 1:
            return "Parabolic"
        elif self.e > 1:
            return "Hyperbolic"
        else:
            return "Invalid"

    def turn_angle(self):
        """
        Calculates the hyperbolic turn angle.
        """
        if self.e <= 1:
            raise ValueError("Can't calculate hyperbolic turn angle for elliptical or circular orbits.")
        return np.degrees(2 * np.arcsin(1/self.e))
    
    def aiming_radius(self):
        """
        Calculates the aiming radius for hyperbolic orbits.
        """
        if self.e <= 1:
            raise ValueError("Can't calculate aiming radius for elliptical or circular orbits.")
        return self.b

    def C3(self):
        """
        Calculates the hyperbolic energy.
        """
        return self.epsilon * 2

    def gamma_at_theta(self, theta):
        """
        Calculates the flight path angle at a given true anomaly.
        """
        theta = np.radians(theta)
        return np.degrees(np.arctan2(self.e * np.sin(theta), 1 + self.e * np.cos(theta)))

    def gamma_at_state_vectors(self, r_vec, v_vec):
        """
        Calcula o ângulo de trajetória (flight path angle) a partir dos vetores de estado.
        
        :param r_vec: Vetor posição (km)
        :param v_vec: Vetor velocidade (km/s)
        :return: Ângulo de trajetória em graus
        """
        # Calcular os módulos dos vetores
        r = np.linalg.norm(r_vec)
        v = np.linalg.norm(v_vec)
        
        # Calcular o produto escalar
        r_dot_v = np.dot(r_vec, v_vec)
        
        # Calcular o ângulo usando arccos do produto escalar normalizado
        gamma = np.arccos(r_dot_v/(r*v))
        
        # Converter para graus e subtrair de 90° para obter o flight path angle
        return 90 - np.degrees(gamma)

    #######   Alternative variables   #######
    def M_from_E(self, E):
        """
        Calculates the mean anomaly using the eccentric anomaly.
        """
        E = np.radians(E)
        if self.e < 1:
            return np.degrees(E - self.e * np.sin(E))
        elif self.e > 1:
            return np.degrees(self.e * np.sinh(E) - E)
        else:
            raise ValueError("Can't calculate mean anomaly from eccentric anomaly for parabolic orbits.")

    def M_from_theta(self, theta):
        """
        Calculates the mean anomaly using the true anomaly.
        """
        if self.e > 1 or self.e < 1:
            raise ValueError("Can't calculate mean anomaly from true anomaly for elliptical or hyperbolic orbits.")
        
        theta = np.radians(theta)
        return np.degrees(1/6 * np.tan(theta/2)**3 + 1/2 * np.tan(theta/2))

    def M_at_t(self, t):
        """
        Calculates the mean anomaly using the time.
        """
        if self.e < 1:
            return np.degrees(2 * np.pi * t / self.T)
        elif self.e == 1:
            return np.degrees(self.mu**2 * t/ self.h**3)
        else:
            return np.degrees(self.mu**2 / self.h**3 * (self.e**2 - 1)**(3/2) * t)

    def E_from_M(self, M):
        """
        Calculates the eccentric anomaly using the mean anomaly.
        """
        if self.e == 1:
            raise ValueError("Can't calculate eccentric anomaly for parabolic orbits.")
        
        M = np.radians(M)
        if self.e < 1:
            if M < np.pi:
                E0 = M + self.e/2
            else:
                E0 = M - self.e/2
        else:
            E0 = M

        def f(E, M):
            if self.e < 1:
                return E - self.e * np.sin(E) - M
            else:
                return self.e * np.sinh(E) - E - M

        def f_prime(E, M):
            if self.e < 1:
                return 1 - self.e * np.cos(E)
            else:
                return self.e * np.cosh(E) - 1
        
        E = root_scalar(f, fprime=f_prime, x0=E0, args=(M), method='newton').root

        return np.degrees(E)

    def E_from_theta(self, theta):
        """
        Calculates the eccentric anomaly using the true anomaly.
        """
        theta = np.radians(theta)
        if self.e < 1:
            return np.degrees(2 * np.arctan(np.sqrt((1 - self.e) / (1 + self.e)) * np.tan(theta/2)))
        elif self.e > 1:
            if np.abs(theta) > np.radians(self.theta_inf()):
                raise ValueError(f"Can't calculate eccentric anomaly for true anomaly greater than hyperbolic infinity angle: {self.theta_inf()} degrees.")
            return np.degrees(2 * np.arctanh(np.sqrt((self.e - 1) / (self.e + 1)) * np.tan(theta/2)))
        else:
            raise ValueError("Can't calculate eccentric anomaly from true anomaly for parabolic orbits.")

    def Qui_at_delta_t(self, delta_t):
        """
        Finds the universal anomaly at a given time.
        """
        Q0 = np.sqrt(self.mu) * np.abs(self.alpha) * delta_t
        
        def f(Qui):
            z = self.alpha * Qui**2
            Qui = (self.r0*self.vr0 / np.sqrt(self.mu) * Qui**2 * self.C(z) + 
                    (1 - self.alpha*self.r0) * Qui**3 * self.S(z) + 
                    self.r0 * Qui - np.sqrt(self.mu) * delta_t)
            return Qui
        def f_dot(Qui):
            z = self.alpha * Qui**2
            Qui_dot = (self.r0*self.vr0 / np.sqrt(self.mu) * Qui * (1 - self.alpha *Qui**2*self.S(z)) +
                    (1 - self.alpha*self.r0) * Qui**2 * self.C(z) + self.r0)
            return Qui_dot
        
        Qui = root_scalar(f, fprime=f_dot, x0=Q0, method='newton').root
        return Qui

    #######   Lagrange coefficients   #######
    def lagrange_coefficients_true(self, r_vec_0, v_vec_0, delta_theta):
        """
        Calculates the Lagrange coefficients using true anomaly.

        :param r_vec_0: Initial position vector (km)
        :param v_vec_0: Initial velocity vector (km/s)
        :param delta_theta: Change in true anomaly (degrees)
        :return: f, g, f_dot, g_dot
        """
        delta_theta = np.radians(delta_theta)
        r0 = np.linalg.norm(r_vec_0)
        v0 = np.linalg.norm(v_vec_0)
        h = np.linalg.norm(np.cross(r_vec_0, v_vec_0))
        vr0 = np.dot(r_vec_0, v_vec_0) / r0

        r = h**2 / (self.mu * (1 + (h**2/(self.mu*r0) - 1)*np.cos(delta_theta) - h*vr0/self.mu * np.sin(delta_theta)))
        f = 1 - self.mu*r/h**2 * (1 - np.cos(delta_theta))
        g = r*r0/h * np.sin(delta_theta)
        f_dot = self.mu/h * (1 - np.cos(delta_theta))/np.sin(delta_theta) * (self.mu/h**2 * (1 - np.cos(delta_theta)) - 1/r0 - 1/r)
        g_dot = 1 - self.mu*r0/h**2 * (1 - np.cos(delta_theta))

        return f, g, f_dot, g_dot

    def lagrange_coefficients_universal(self, r_vec_0, v_vec_0, Qui, delta_t):
        """
        Calculates the Lagrange coefficients using universal anomaly.

        :param r_vec_0: Initial position vector (km)
        :param v_vec_0: Initial velocity vector (km/s)
        :param Qui: Universal anomaly (km^1/2)
        :param delta_t: Time variation (s)
        :return: f, g, f_dot, g_dot
        """
        r0 = np.linalg.norm(r_vec_0)
        v0 = np.linalg.norm(v_vec_0)
        vr0 = np.dot(r_vec_0, v_vec_0) / r0
        alpha = self.alpha_from_mu_r_v(self.mu, r0, v0)
        z = alpha * Qui**2

        f = 1 - Qui**2/r0 * self.C(z)
        g = delta_t - Qui**3/np.sqrt(self.mu) * self.S(z)
        
        r_vec = f*r_vec_0 + g*v_vec_0
        r = np.linalg.norm(r_vec)

        f_dot = np.sqrt(self.mu)/(r*r0) * (alpha*Qui**3 * self.S(z) - Qui)
        g_dot = 1 - Qui**2/r * self.C(z)

        return f, g, f_dot, g_dot


    #######   Distance   #######
    def r_avg_by_theta(self):
        """
        Calculates the theta-averaged radius using periapsis and apoapsis distances.
        Note: Can't be used for parabolic or hyperbolic orbits.
        """
        if self.e >= 1:
            raise ValueError("Can't calculate theta-averaged radius for parabolic or hyperbolic orbits.")
        return (self.rp * self.ra)**(1/2)

    def r_avg_by_t(self):
        """
        Calculates the time-averaged radius using semi-major axis and eccentricity.
        Note: Can't be used for parabolic or hyperbolic orbits.
        """
        if self.e >= 1:
            raise ValueError("Can't calculate time-averaged radius for parabolic or hyperbolic orbits.")
        return self.a * (1 + self.e**2/2)

    def r_at_theta(self, theta):
        """
        Calculates the distance from the primary body center to a point on the orbit at a given true anomaly.
        """
        theta = np.radians(theta)
        return self.h**2 / (self.mu) * 1/(1 + self.e * np.cos(theta))


    #######   Velocities   #######
    def v_esc(self, r):
        """
        Calculates the escape velocity at a given distance.
        """
        return np.sqrt(2 * self.mu / r)

    def v_inf(self):
        """
        Calculates the hyperbolic excess velocity.
        """
        if self.e <= 1:
            raise ValueError("Can't calculate hyperbolic excess velocity for elliptical or parabolic orbits.")
        return np.sqrt(self.mu/self.a)

    def v_at_r(self, r):
        """
        Calculates the orbital velocity at a given distance using the vis-viva equation.
        
        :param r: Distance from the primary body center to the point (km)
        :return: Velocity at the point (km/s)
        """
        if self.e == 1:
            return np.sqrt(self.mu * 2 / r)
        else:
            return np.sqrt(self.mu * (2 / r - 1 / self.a))


    #######   Theta evolution   #######
    def theta_inf(self):
        """
        Calculates the true anomaly at infinity for hyperbolic orbits.
        """
        if self.e <= 1:
            raise ValueError("Can't calculate true anomaly at infinity for elliptical or circular orbits.")
        return np.degrees(np.arccos(-1/self.e))

    def theta_from_r_gamma(self, r, gamma):
        """
        Calculates the true anomaly using the flight path angle.
        """
        gamma = np.radians(gamma)
        if self.e == 0:  # circular orbit
            theta = 0 if gamma >= 0 else 180
        else:
            sin_theta = np.sin(gamma) * (1 + self.h**2/(self.mu*r)) / self.e
            cos_theta = (self.h**2/(self.mu*r) - 1) / self.e
            theta = np.mod(np.degrees(np.arctan2(sin_theta, cos_theta)), 360)
        return theta

    def theta_from_E(self, E):
        """
        Calculates the true anomaly using the eccentric anomaly.
        """
        E = np.radians(E)
        if self.e < 1:
            return np.degrees(2 * np.arctan(np.sqrt((1 + self.e) / (1 - self.e)) * np.tan(E/2)))
        elif self.e > 1:
            return np.degrees(2 * np.arctan(np.sqrt((self.e + 1) / (self.e - 1)) * np.tanh(E/2)))
        else:
            raise ValueError("Can't calculate true anomaly from eccentric anomaly for parabolic orbits.")

    def theta_at_r(self, r):
        """
        Calculates the true anomaly at a given distance from the primary body center.
        Returns both possible values as a list, with angles between 0 and 360 degrees.
        """
        theta1 = np.arccos((self.h**2 / (self.mu * r) - 1) / self.e)
        theta2 = 2*np.pi - theta1
        return [np.degrees(theta1), np.degrees(theta2)]

    def theta_at_t(self, t):
        """
        Calculates the true anomaly at a given time.
        """
        if self.e < 1 or self.e > 1: # elliptical or hyperbolic
            M = self.M_at_t(t)
            E = self.E_from_M(M)
            return np.mod(self.theta_from_E(E), 360)
        else: # parabolic
            M = self.M_at_t(t)
            M = np.radians(M)
            z = (3*M + np.sqrt(1 + 9*M**2))**(1/3)
            theta = np.degrees(2 * np.arctan(z - 1/z))
            return np.mod(theta, 360)

    def theta_at_state_vectors(self, r_vec, v_vec):
        """
        Calculates the true anomaly at a given state vectors.
        """
        h_vec = np.cross(r_vec, v_vec)
        e_vec = np.cross(v_vec, h_vec)/self.mu - r_vec/np.linalg.norm(r_vec)
        if self.e == 0:
            theta = 0
        else:
            cos_theta = np.dot(e_vec, r_vec)/(self.e * np.linalg.norm(r_vec))
            sin_theta = np.dot(h_vec, np.cross(e_vec, r_vec))/(self.e * self.h * np.linalg.norm(r_vec))
            theta = np.mod(np.degrees(np.arctan2(sin_theta, cos_theta)), 360)
        return theta


    #######   Time evolution   #######
    def t_at_M(self, M):
        """
        Calculates the time using the mean anomaly.
        """
        M = np.radians(M)
        if self.e < 1:
            t = M * self.T / (2 * np.pi)
            while t < 0:
                t += self.T
        elif self.e == 1:
            t = M * self.h**3 / self.mu**2
        else:
            t = M * self.h**3 / (self.mu**2 * (self.e**2 - 1)**(3/2))
        return t
    
    def t_at_theta(self, theta):
        """
        Calculates the time using the true anomaly.
        """
        if self.e < 1 or self.e > 1: # elliptical or hyperbolic
            E = self.E_from_theta(theta)
            M = self.M_from_E(E)
            return self.t_at_M(M)
        else: # parabolic
            M = self.M_from_theta(theta)
            return self.t_at_M(M)


    #######   State vectors   #######
    def r_vec_at_theta(self, theta, frame="perifocal"):
        """
        Calculates the position vector at a given true anomaly, in a given frame of reference.
        """
        theta_rad = np.radians(theta)
        if frame == "perifocal":
            return self.h**2 / (self.mu * (1 + self.e*np.cos(theta_rad))) * np.array([np.cos(theta_rad), np.sin(theta_rad), 0])
        elif frame == "bodycentric":
            return self.Q_bc_peri.T @ self.r_vec_at_theta(theta, frame="perifocal")
        elif frame == "rotatingBodycentric":
            return OrbitalBase.convert_bc_to_rbc(self.r_vec_at_theta(theta, frame="bodycentric"), self.omega_body)
        else:
            raise ValueError("Invalid frame of reference.")
    
    def v_vec_at_theta(self, theta, frame="perifocal"):
        """
        Calculates the velocity vector at a given true anomaly, in a given frame of reference.
        """
        theta_rad = np.radians(theta)
        if frame == "perifocal":
            return self.mu / self.h * np.array([-np.sin(theta_rad), self.e + np.cos(theta_rad), 0])
        elif frame == "bodycentric":
            return self.Q_bc_peri.T @ self.v_vec_at_theta(theta, frame="perifocal")
        elif frame == "rotatingBodycentric":
            return self.convert_bc_to_rbc(self.v_vec_at_theta(theta, frame="bodycentric"), self.omega_body)
        else:
            raise ValueError("Invalid frame of reference.")    

    def state_vectors_at_Qui(self, Qui, delta_t):
        """
        Finds the position and velocity vectors at a given universal anomaly and delta time.
        """
        z = self.alpha * Qui**2
        
        ## Lagrange coefficients ##
        f = 1 - Qui**2/self.r0 * self.C(z)
        g = delta_t - Qui**3/np.sqrt(self.mu) * self.S(z)
        
        r_vec = f*self.r_vec_0 + g*self.v_vec_0
        r = np.linalg.norm(r_vec)

        ## Lagrange coefficients derivatives ##
        f_dot = np.sqrt(self.mu)/(r*self.r0) * (self.alpha*Qui**3 * self.S(z) - Qui)
        g_dot = 1 - Qui**2/r * self.C(z)

        v_vec = f_dot*self.r_vec_0 + g_dot*self.v_vec_0

        return r_vec, v_vec

    def state_vectors_at_theta(self, theta, frame="perifocal"):
        """
        Calculates the state vector at a given true anomaly, in a given frame of reference.
        """
        if frame == "perifocal":
            r_vec = self.r_vec_at_theta(theta, frame)
            v_vec = self.v_vec_at_theta(theta, frame)
        elif frame == "bodycentric":
            r_vec = self.Q_bc_peri.T @ self.r_vec_at_theta(theta, frame="perifocal")
            v_vec = self.Q_bc_peri.T @ self.v_vec_at_theta(theta, frame="perifocal")
        elif frame == "rotatingBodycentric":
            r_vec = self.r_vec_at_theta(theta, frame="rotatingBodycentric")
            v_vec = self.v_vec_at_theta(theta, frame="rotatingBodycentric")
        else:
            raise ValueError("Invalid frame of reference.")
        return r_vec, v_vec

    def state_vectors_at_t(self, t):
        """
        Finds the position and velocity vectors at a given time.
        """
        if self.frame == "perifocal":
            delta_t = t - self.t0_clock
            Qui = self.Qui_at_delta_t(delta_t)
            r_vec, v_vec = self.state_vectors_at_Qui(Qui, delta_t)
        elif self.frame == "bodycentric":
            delta_t = t - self.t0_clock
            Qui = self.Qui_at_delta_t(delta_t)
            r_vec, v_vec = self.state_vectors_at_Qui(Qui, delta_t)
            r_vec = self.convert_bodycentric_to_rotatingBodycentric(r_vec, self.omega_body)
            v_vec = self.convert_bodycentric_to_rotatingBodycentric(v_vec, self.omega_body)
        return r_vec, v_vec

    def state_after_delta_theta(self, r_vec_0, v_vec_0, delta_theta):
        """
        Calculates the position and velocity vectors after a given true anomaly change.

        :param r_vec_0: Initial position vector (km)
        :param v_vec_0: Initial velocity vector (km/s)
        :param delta_theta: Change in true anomaly (degrees)
        :return: Position and velocity vectors after the change (km and km/s)
        """
        f, g, f_dot, g_dot = self.lagrange_coefficients_true(r_vec_0, v_vec_0, delta_theta)
        r_vec = f*r_vec_0 + g*v_vec_0
        v_vec = f_dot*r_vec_0 + g_dot*v_vec_0
        return r_vec, v_vec
    
    def state_after_delta_t(self, r_vec_0, v_vec_0, delta_t):
        """
        Calculates the position and velocity vectors after a given time change.
        """
        Qui = self.Qui_at_delta_t(delta_t)
        return self.state_vectors_at_Qui(Qui, delta_t)

    #######   Adding Points   #######
    def add_orbital_position(self, theta, t=None, name=None):
        """
        Adds an orbital position to the orbit at a given true anomaly and (optionally) time.

        :param theta: True anomaly (degrees)
        :param t: Time (s), default is None
        :param name: Name of the position, default is None
        """
        r = self.r_at_theta(theta)
        v = self.v_at_r(r)
        position = {
            'name': name,
            'r': r,
            'v': v, 
            'theta': theta,
            't': t
        }
        
        self.positions.append(position)
    
    def add_point_in_orbital_plane(self, r, theta, name=None):
        """
        Adds a point in the orbital plane at a given distance and true anomaly.
        """
        point = {
            'name': name,
            'r': r,
            'theta': theta
        }
        
        self.points_in_orbital_plane.append(point)

    def add_zero_state_from_theta(self, theta0=0, t0_clock=0, frame="perifocal", add_position=True):
        """
        Adds the initial state vectors for trajectory calculations.
        """
        self.theta0 = theta0
        r_vec_0, v_vec_0 = self.state_vectors_at_theta(theta0, frame=frame)
        self.r_vec_0 = r_vec_0
        self.v_vec_0 = v_vec_0
        self.r0 = np.linalg.norm(r_vec_0)
        self.v0 = np.linalg.norm(v_vec_0)
        self.vr0 = np.dot(r_vec_0, v_vec_0)/self.r0
        self.t0_clock = t0_clock
        self.t0_orbit = self.t_at_theta(self.theta0)
        self.delta_t = self.t0_clock - self.t0_orbit
        if add_position:
            self.add_orbital_position(self.theta0, self.t0_clock, name="Initial Position")

    def add_zero_state_from_vectors(self, r_vec_0, v_vec_0, t0_clock=0, add_position=True):
        """
        Adds the initial state vectors for trajectory calculations.
        """
        self.r_vec_0 = r_vec_0
        self.v_vec_0 = v_vec_0
        self.r0 = np.linalg.norm(r_vec_0)
        self.v0 = np.linalg.norm(v_vec_0)
        self.vr0 = np.dot(r_vec_0, v_vec_0)/self.r0
        self.theta0 = self.theta_at_state_vectors(r_vec_0, v_vec_0)
        self.t0_clock = t0_clock
        self.t0_orbit = self.t_at_theta(self.theta0)
        self.delta_t = self.t0_clock - self.t0_orbit
        if add_position:
            self.add_orbital_position(self.theta0, self.t0_clock, name="Initial Position")

    def add_trajectory_points(self, r, theta, t_clock=None):
        """
        Adds a trajectory in the orbital plane between two given true anomalies.
        """

        point = {
            'r': r,
            'theta': theta,
            't_clock': t_clock
        }
        
        self.trajectory_points.append(point)
    
    def trajectory(self, theta1=360, t1_clock=None, n_points=20, add_trajectory_points=True):
        """
        Calculates the trajectory between two true anomalies or times.

        :param theta_1: Final true anomaly (degrees) (optional, default is 360)
        :param t1: Final time (s) - if provided, overrides theta_1 (optional, default is None)
        :param n_points: Number of points to evaluate the trajectory (optional, default is 20)
        :param add_trajectory_points: Boolean to add the trajectory points to the orbit (optional, default is True)
        """
        
        # Minimum interval (s) between points
        min_interval = 10
        
        if t1_clock is None:
            # If t1 is not provided, calculate from theta1
            if self.e < 1: # elliptical and circular
                # Calculate number of complete orbits
                n_orbits = int(theta1 / 360)
                
                # Calculate angle within first orbit
                theta1_norm = theta1 % 360
                
                # Calculate time for normalized angle
                self.t1_orbit = self.t_at_theta(theta1_norm)
                
                # Add time for complete orbits
                self.t1_orbit += n_orbits * self.T
                
                # Ensure t1_orbit is greater than t0_orbit
                while self.t1_orbit < self.t0_orbit:
                    self.t1_orbit += self.T
            else: # parabolic and hyperbolic
                self.t1_orbit = self.t_at_theta(theta1)
            
            self.t1_clock = self.delta_t + self.t1_orbit
        else:
            self.t1_clock = t1_clock
            self.t1_orbit = self.t1_clock - self.delta_t
        
        

        # Add final position
        theta1 = self.theta_at_t(self.t1_orbit)
        self.add_orbital_position(theta1, self.t1_clock, name='Final Position')
        
        # Ensure interval between points is not smaller than min_interval
        n_points = min(n_points, max(2, int((self.t1_clock - self.t0_clock) / min_interval)))
        t_eval = np.linspace(self.t0_clock, self.t1_clock, n_points)

        # Calculate trajectory for each evaluation point
        theta_prev = self.theta0
        for t_clock in t_eval:
            t_orbit = t_clock - self.delta_t
            theta = self.theta_at_t(t_orbit)
            if theta < theta_prev:
                theta += 360
            r = self.r_at_theta(theta)
            if add_trajectory_points:
                self.add_trajectory_points(r, theta, t_clock)
            theta_prev = theta

    #######   Plotting   #######
    def plot(self, orbit=True, points=True, positions=True, velocities=True, trajectory=True, frame="perifocal", plot3d=True):
        """
        Plots the orbit with optional points and velocities.

        :param orbit: Boolean to plot the orbit
        :param points: Boolean to plot the points
        :param positions: Boolean to plot the already added orbital positions
        :param velocities: Boolean to plot the velocities
        :param trajectory: Boolean to plot the trajectory
        :param frame: Frame of reference to plot the orbit
        """
        fig = plt.figure(figsize=(16, 9))
        if plot3d:
            ax = fig.add_subplot(projection='3d')
        else:
            ax = fig.add_subplot()
        # List of colors
        colors = ['red', 'orange', 'purple', 'brown', 'deeppink', 'olive', 'deepskyblue', 'blue', 'green']

        def draw_orbit():
            scale = self.rp  # Scale for auxiliary vectors
            # Angular momentum vector
            h_vec = self.h_vec/self.h * scale #frame bodycentric
            # Eccentricity vector
            e_vec = self.e_vec/self.e * scale #frame bodycentric
            if frame == "perifocal":
                h_vec = self.Q_bc_peri @ h_vec
                e_vec = self.Q_bc_peri @ e_vec
            
            # Create array of true anomaly
            if self.e == 1:
                theta_min = -120
                theta_max = 120
                if positions:
                    for pos in self.positions:
                        if pos['theta'] < theta_min:
                            theta_min = pos['theta']
                        if pos['theta'] > theta_max:
                            theta_max = pos['theta']
                theta = np.linspace(theta_min, theta_max, 360)
            elif self.e > 1:
                epsilon = 15
                theta_min = -self.theta_inf() + epsilon
                theta_max = self.theta_inf() - epsilon
                if positions:
                    for pos in self.positions:
                        if pos['theta'] < theta_min:
                            theta_min = pos['theta']
                        if pos['theta'] > theta_max:
                            theta_max = pos['theta']
                theta = np.linspace(theta_min, theta_max, 360)
            else:
                theta = np.linspace(0, 360, 360)
            
            r_vecs = np.array([self.r_vec_at_theta(t, frame) for t in theta])
            
            # Plot orbit
            if plot3d:
                ax.plot(r_vecs[:,0], r_vecs[:,1], r_vecs[:,2], 'b-', label='Orbit', linewidth=2)
                ax.quiver(0, 0, 0, h_vec[0], h_vec[1], h_vec[2], color='g', label='h')
                ax.quiver(0, 0, 0, e_vec[0], e_vec[1], e_vec[2], color='r', label='e')
            else:
                ax.plot(r_vecs[:,0], r_vecs[:,1], 'b-', label='Orbit', linewidth=2, zorder=1)
                ax.quiver(0, 0, h_vec[0], h_vec[1], color='g', label='h', scale_units='xy', scale=1, width=0.003)
                ax.quiver(0, 0, e_vec[0], e_vec[1], color='r', label='e', scale_units='xy', scale=1, width=0.003)

        def draw_body(x_center, y_center, z_center=0, index=1, radius=None, color=None, label=None):
            if index == 1:
                if radius is None:
                    radius = self.body1radius
                if color is None:
                    color = 'red'
                if label is None:
                    label = 'Central Body'
            elif index == 2:
                if radius is None:
                    radius = self.body2radius
                if color is None:
                    color = 'blue'
                if label is None:
                    label = 'Secondary Body'
            
            if plot3d:
                if radius is not None:
                    u = np.linspace(0, 2 * np.pi, 36)
                    v = np.linspace(0, np.pi, 36)
                    x = x_center + radius * np.outer(np.cos(u), np.sin(v))
                    y = y_center + radius * np.outer(np.sin(u), np.sin(v))
                    z = z_center + radius * np.outer(np.ones(np.size(u)), np.cos(v))
                    ax.plot_surface(x, y, z, color=color, alpha=0.4, label=label)
                else:
                    ax.plot(x_center, y_center, z_center, 'o', color=color, label=label)
            else:
                if radius is not None:
                    circle = plt.Circle((x_center, y_center), radius, color=color, alpha=0.4, label=label, zorder=2)
                    plt.gca().add_patch(circle)
                else:
                    plt.plot(x_center, y_center, 'o', color=color, label=label, zorder=2)
            
        def draw_points():
            if self.points_in_orbital_plane:
                for i, point in enumerate(self.points_in_orbital_plane):
                    label = point['name']
                    color = 'gray'
                    r_vec = self.r_vec_at_theta(point['theta'], frame=frame)
                    if plot3d:
                        ax.plot(r_vec[0], r_vec[1], r_vec[2], 'o', color=color, label=label)
                    else:
                        ax.plot(r_vec[0], r_vec[1], 'o', color=color, label=label, zorder=2)
        
        def draw_positions():
            if self.positions:
                for i, position in enumerate(self.positions):
                    label = position.get('name', f'Position {i+1}')
                    color = colors[i % len(colors)]  # Recycle colors if there are more points than colors
                    r_vec = self.r_vec_at_theta(position['theta'], frame=frame)
                    if position['t'] is not None:
                        label += f' (t={position["t"]:.2f}s)'
                    # Plot the secondary body if provided
                    draw_body(r_vec[0], r_vec[1], r_vec[2], index=2, radius=self.body2radius, color=color, label=label)
                    
                    # Plot velocity vector
                    if velocities:
                        v_vec = self.v_vec_at_theta(position['theta'], frame=frame)
                        v_vec = v_vec * self.p/20  # Velocity vector scale
                        if plot3d:
                            ax.quiver(r_vec[0], r_vec[1], r_vec[2], v_vec[0], v_vec[1], v_vec[2], color='r')
                        else:
                            ax.quiver(r_vec[0], r_vec[1], v_vec[0], v_vec[1], color='r', scale_units='xy', scale=1, width=0.003, zorder=5)

        def draw_trajectory():
            if self.trajectory_points:
                x_points = []
                y_points = []
                z_points = []
                annotations = []  # List to store all annotations

                # Trajectory points for annotations
                for point in self.trajectory_points:
                    r_vec = self.r_vec_at_theta(point['theta'], frame=frame)
                    x_points.append(r_vec[0])
                    y_points.append(r_vec[1])
                    z_points.append(r_vec[2])
                    if point['t_clock'] is not None:
                        t_min = min(p['t_clock'] for p in self.trajectory_points if p['t_clock'] is not None)
                        t_max = max(p['t_clock'] for p in self.trajectory_points if p['t_clock'] is not None)
                        
                        t_clock = point['t_clock']
                        t_norm = (t_clock - t_min) / (t_max - t_min)
                        
                        r = 1 - t_norm
                        g = 0.5
                        b = t_norm
                        
                        if plot3d:
                            point_plot = ax.plot(r_vec[0], r_vec[1], r_vec[2], '.', color=(r,g,b),
                                                picker=True, pickradius=5)[0]
                        else:
                            point_plot = plt.plot(r_vec[0], r_vec[1], '.', color=(r,g,b),
                                                picker=True, pickradius=5, zorder=4)[0]
                        
                        # Create annotation
                        annotation = plt.annotate(f't = {t_clock:.2f}s',
                            xy=(r_vec[0], r_vec[1]), xytext=(10, 10),
                            textcoords='offset points',
                            bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                            arrowprops=dict(arrowstyle='->'),
                            visible=False,
                            zorder=10)
                        
                        annotations.append((point_plot, annotation))
                
                # Hover function that manages all annotations
                def hover(event):
                    if event.inaxes == plt.gca():
                        for point_plot, annotation in annotations:
                            cont, _ = point_plot.contains(event)
                            annotation.set_visible(cont)
                        plt.draw()
                
                plt.gcf().canvas.mpl_connect('motion_notify_event', hover)
                
                # Generate high-resolution points for the continuous line
                theta_min = min(p['theta'] for p in self.trajectory_points)
                theta_max = max(p['theta'] for p in self.trajectory_points)
                thetas = np.linspace(theta_min, theta_max, int(theta_max - theta_min))
                x_line = []
                y_line = []
                z_line = []
                for theta in thetas:
                    r_vec = self.r_vec_at_theta(theta, frame=frame)
                    x_line.append(r_vec[0])
                    y_line.append(r_vec[1])
                    z_line.append(r_vec[2])
                
                if plot3d:
                    ax.plot(x_line, y_line, z_line, '-', color='black', label='Trajectory')
                else:
                    plt.plot(x_line, y_line, '-', color='black', label='Trajectory', zorder=3)

        if frame == "perifocal" or frame == "bodycentric":
            ### Plot main body ###
            draw_body(0, 0, 0, index=1)

            ### Plot orbit ###
            if orbit:
                draw_orbit()
            
            # Plot points
            if points:
                draw_points()

            ### Plot positions ###
            if positions:
                draw_positions()
            
            ### Plot trajectory points ###
            if trajectory:
                draw_trajectory()

        elif frame == "rotatingBarycentric":
            plt.title('Frame: Rotating Barycentric')

            ### Plot main bodies ###
            draw_body(self.body1_center[0], self.body1_center[1], index=1)
            draw_body(self.body2_center[0], self.body2_center[1], index=2)
            
            ### Plot points ###
            if points:
                draw_points()

            ### Plot positions ###
            if positions:
                draw_positions()
            
            ### Plot trajectory points ###
            if trajectory:
                draw_trajectory()
        
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        if plot3d:
            ax.set_zlabel('Z (km)')
        ax.axis('equal')
        ax.grid(True)
        plt.title(f'Frame: {frame}')
        plt.tight_layout()
        plt.get_current_fig_manager().window.showMaximized()
        plt.show()

    ####### String representation #######
    def __str__(self):
        """
        Returns a string representation of the orbital parameters.
        """
        # Get all attributes of the class
        atributos = vars(self)
        # Initial string
        resultado = "Orbital Parameters:\n"
        # Add each attribute to the string
        for nome, valor in atributos.items():
            # Skip attributes that are None or that should not be shown
            if valor is None or nome in ['points_in_orbital_plane', 'positions', 'trajectory_points']:
                continue
                
            # Format numbers as float with 3 decimal places
            if isinstance(valor, (int, float)):
                resultado += f"{nome}: {valor:.3f}\n"
            else:
                resultado += f"{nome}: {valor}\n"

        return resultado




class Orbit(OrbitalBase):
    def __init__(self, mu=None, m1=None, m2=0, a=None, e=None, rp=None, ra=None, h=None, Omega=None, i=None, omega=None, theta=None, t0_clock=0, body1radius=None, body2radius=None, **orbital_params):
        """
        Initializes the Orbit class with primary orbital parameters.
        :param mu: Standard gravitational parameter (km³/s²)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param a: Semi-major axis (km)
        :param e: Eccentricity
        :param rp: Pericenter radius (km)
        :param ra: Apocenter radius (km)
        :param h: Angular momentum (km²/s)
        :param body1radius: Radius of the central body (km) (optional)
        :param body2radius: Radius of the orbiting body (km) (optional)
        :param orbital_params: Optional dictionary with pre-calculated orbital parameters (when initializing from state vectors)
        """
        super().__init__()
        self.body1radius = body1radius
        self.body2radius = body2radius

        # If orbital_params is provided, first add explicit parameters to the dictionary
        if orbital_params:
            # Add explicit parameters that aren't None to orbital_params
            explicit_params = {
                'mu': mu,
                'm1': m1,
                'm2': m2,
                'a': a,
                'e': e,
                'rp': rp,
                'ra': ra,
                'h': h,
                'Omega': Omega,
                'i': i,
                'omega': omega,
                'theta': theta,
                't0_clock': t0_clock
            }
            # Update orbital_params only with non-None values
            orbital_params.update({k: v for k, v in explicit_params.items() if v is not None})
            
            # Now set all attributes
            for key, value in orbital_params.items():
                setattr(self, key, value)
            return
        else:
            # Calculate mu if not provided
            if mu is not None:
                self.mu = mu  # Gravitational parameter
            elif m1 is not None:
                self.mu = self.mu_from_m1_m2(m1, m2)
            else:
                raise ValueError("Provide either mu or m1 (and optionally m2) to calculate mu.")
            self.m1 = m1
            self.m2 = m2
            
            # Calculate rp, e, h, ra, a, b, epsilon, T, p, alpha
            if rp is not None and e is not None:
                self.rp = rp
                self.e = e
                self.h = self.h_from_mu_rp_e(self.mu, self.rp, self.e)
                self.ra = self.ra_from_mu_h_e(self.mu, self.h, self.e)
                self.a = self.a_from_rp_ra(self.rp, self.ra)
            elif rp is not None and ra is not None:
                self.rp = rp
                self.ra = ra
                self.a = self.a_from_rp_ra(self.rp, self.ra)
                self.e = self.e_from_rp_ra(self.rp, self.ra)
                self.h = self.h_from_mu_rp_e(self.mu, self.rp, self.e)
            elif rp is not None and h is not None:
                self.rp = rp
                self.h = h
                self.e = self.e_from_mu_h_rp(self.mu, self.h, self.e)
                self.ra = self.ra_from_mu_h_e(self.mu, self.h, self.e)
                self.a = self.a_from_rp_ra(self.rp, self.ra)
            elif e is not None and h is not None:
                self.e = e
                self.h = h
                self.rp = self.rp_from_mu_h_e(self.mu, self.h, self.e)
                self.ra = self.ra_from_mu_h_e(self.mu, self.h, self.e)
                self.a = self.a_from_rp_ra(self.rp, self.ra)
            elif a is not None and e is not None:
                if e == 1:
                    raise ValueError("Can't find rp from 'a' and 'e' for parabolic (e=1) orbits.")
                self.a = a
                self.e = e
                self.rp = self.rp_from_a_e(self.a, self.e)
                self.ra = self.ra_from_a_e(self.a, self.e)
                self.h = self.h_from_mu_rp_e(self.mu, self.rp, self.e)
            else:
                raise ValueError("Provide sufficient parameters to define the orbit (rp and e, rp and ra, rp and h, e and h or a and e).")
            
            self.b = self.b_from_a_e(self.a, self.e)
            self.epsilon = self.epsilon_from_mu_h_e(self.mu, self.h, self.e)
            self.T = self.T_from_mu_a(self.mu, self.a)
            self.p = self.p_from_mu_h(self.mu, self.h)
            self.alpha = self.alpha_from_a(self.a)

            if Omega is not None and i is not None and omega is not None and theta is not None:
                self.Omega = Omega
                self.i = i
                self.omega = omega
                self.theta = theta
                self.Q_bc_peri = self.Q_from_Euler_angles(Omega, i, omega, pattern="classical")
                self.r_vec, self.v_vec = self.state_vectors_at_theta(theta, frame="bodycentric")
                self.h_vec = self.h_vec_from_mu_state_vectors(self.mu, self.r_vec, self.v_vec)
                self.e_vec = self.e_vec_from_mu_state_vectors(self.mu, self.r_vec, self.v_vec)
                self.n_vec = self.n_vec_from_mu_h_vec(self.mu, self.h_vec)
                self.r = np.linalg.norm(self.r_vec)
                self.v = np.linalg.norm(self.v_vec)
                self.n = np.linalg.norm(self.n_vec)

                # Initial parameters
                self.t0_clock = t0_clock
                self.t0_orbit = self.t_at_theta(self.theta)
                self.delta_t = self.t0_clock - self.t0_orbit
                self.r_vec_0 = self.r_vec
                self.v_vec_0 = self.v_vec
                self.r0 = np.linalg.norm(self.r_vec_0)
                self.v0 = np.linalg.norm(self.v_vec_0)
                self.vr0 = np.dot(self.r_vec_0, self.v_vec_0)/self.r0
                self.theta0 = self.theta_at_state_vectors(self.r_vec_0, self.v_vec_0)
                self.add_orbital_position(self.theta0, self.t0_clock, name="Initial Position")

    @classmethod
    def init_from_2positions(cls, r1, theta1, r2, theta2, mu=None, m1=None, m2=0, body1radius=None, body2radius=None):
        """
        Creates a new instance of Orbit from two orbital positions (r,θ).
        
        :param r1: Radius of the first position (km)
        :param theta1: True anomaly of the first position (degrees)
        :param r2: Radius of the second position (km)
        :param theta2: True anomaly of the second position (degrees)
        :param mu: Standard gravitational parameter (km³/s²)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param body1radius: Radius of the central body (km) (optional)
        :param body2radius: Radius of the orbiting body (km) (optional)
        :return: New instance of Orbit
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Provide mu or m1 (and optionally m2) to calculate mu.")
        
        theta1_rad = np.radians(theta1)
        theta2_rad = np.radians(theta2)
        
        e = (r2/r1 - 1) / (np.cos(theta1_rad) - (r2/r1)*np.cos(theta2_rad))
        h = np.sqrt(mu * r1 * (1 + e * np.cos(theta1_rad)))
        
        orbit_instance = cls(mu=mu, m1=m1, m2=m2, e=e, h=h, body1radius=body1radius, body2radius=body2radius)
        orbit_instance.add_orbital_position(theta=theta1, name="Position 1")
        orbit_instance.add_orbital_position(theta=theta2, name="Position 2")

        return orbit_instance

    @classmethod
    def init_from_r_v_gamma(cls, r, v, gamma, mu=None, m1=None, m2=0, body1radius=None, body2radius=None):
        """
        Creates a new instance of Orbit from distance, velocity and flight path angle.

        :param r: Distance from the primary body center to the position (km)
        :param v: Velocity at the position (km/s)
        :param gamma: Flight path angle (degrees)
        :param mu: Standard gravitational parameter (km³/s²)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param body1radius: Radius of the central body (km) (optional)
        :param body2radius: Radius of the orbiting body (km) (optional)
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Provide mu or m1 (and optionally m2) to calculate mu.")
        
        gamma = np.radians(gamma)
        vr = v * np.sin(gamma)
        vt = v * np.cos(gamma)
        h = r * vt
        e = np.sqrt(vr**2*h**2/mu**2 + (h**2/(mu*r) - 1)**2)
        
        orbit_instance = cls(mu=mu, m1=m1, m2=m2, e=e, h=h, body1radius=body1radius, body2radius=body2radius)
        theta = orbit_instance.theta_from_r_gamma(r, gamma)
        orbit_instance.add_orbital_position(theta=theta, name="Position 1")
        
        return orbit_instance

    @classmethod
    def init_from_state_vectors(cls, r_vec, v_vec, mu=None, m1=None, m2=0, name=None, starting_point=True, t0_clock=0, body1radius=None, body2radius=None):
        """
        Given a position and velocity vectors represented in some inertial frame of reference
        centered in the primary body, creates a new instance of Orbit.

        :param r_vec: Position vector in the inertial frame of reference (km)
        :param v_vec: Velocity vector in the inertial frame of reference (km/s)
        :param mu: Standard gravitational parameter (km³/s²)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param name: Name of the body
        :param starting_point: If True, adds the starting point to the orbit
        :param t0_clock: Time of the starting point (s)
        :param body1radius: Radius of the central body (km) (optional)
        :param body2radius: Radius of the orbiting body (km) (optional)
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Provide mu or m1 (and optionally m2) to calculate mu.")
        
        # Calculate all orbital parameters
        orbital_parameters = cls.orbital_parameters_from_state_vectors(r_vec, v_vec, mu)
        
        # Add parameters that are not in orbital_parameters
        additional_params = {
            'body1radius': body1radius,
            'body2radius': body2radius
        }
        
        # Create the instance passing all the calculated parameters
        orbit_instance = cls(**orbital_parameters, **additional_params)
        
        # Add the starting point if requested
        if starting_point:
            orbit_instance.add_zero_state_from_vectors(r_vec, v_vec, t0_clock)
        else:
            orbit_instance.add_orbital_position(theta=orbital_parameters['theta'], name="Position 1")
        
        return orbit_instance

    @classmethod
    def init_from_r_v_theta(cls, r, v, theta, mu=None, m1=None, m2=0, body1radius=None, body2radius=None):
        """
        Creates a new instance of Orbit from distance, velocity and true anomaly.

        :param r: Distance from the primary body center to the point (km)
        :param v: Velocity at the point (km/s)
        :param theta: True anomaly (degrees)
        :param mu: Standard gravitational    parameter (km³/s²)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param body1radius: Radius of the central body (km) (optional)
        :param body2radius: Radius of the orbiting body (km) (optional)
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Provide mu or m1 (and optionally m2) to calculate mu.")
        
        theta = np.radians(theta)

        e = cls.e_from_mu_r_v_theta(mu, r, v, theta)
        h = cls.h_from_mu_r_e_theta(mu, r, e, theta)

        orbit_instance = cls(mu=mu, m1=m1, m2=m2, e=e, h=h, body1radius=body1radius, body2radius=body2radius)
        orbit_instance.add_orbital_position(theta=theta, name="Position 1")
        
        return orbit_instance




class Three_body_restricted(OrbitalBase):
    """
    Class to calculate the position and velocity vectors of a body in a three-body restricted problem.
    Frame of reference: noninertial, rotating system, with origin at the center of mass and X axis towards m2.
    Note: Valid only for the case where the two main bodies (m1 and m2) are in circular orbit around each other.
    """

    def __init__(self, m1, m2, r12, body1radius=None, body2radius=None):
        """
        Initializes the Three_body_restricted class.

        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg)
        :param r12: Distance between the two bodies (km)
        :param body1radius: Radius of the primary body (km), default is None
        :param body2radius: Radius of the secondary body (km), default is None
        """
        super().__init__()
        self.m1 = m1
        self.m2 = m2
        self.mu1 = self.G * m1
        self.mu2 = self.G * m2
        self.pi1 = m1/(m1 + m2)
        self.pi2 = m2/(m1 + m2)
        self.mu = self.mu1 + self.mu2
        self.r12 = r12
        self.body1_center = np.array([-self.pi2*self.r12, 0, 0])
        self.body2_center = np.array([(1-self.pi2)*self.r12, 0, 0])
        self.body1radius = body1radius
        self.body2radius = body2radius
        self.omega = np.sqrt(self.mu/self.r12**3)

    def lagrange_points(self, add_points=True):
        """
        Calculates the Lagrange points.
        """
        
        #### L4 and L5 ####
        x = self.r12/2 - self.pi2 * self.r12
        y = np.sqrt(3)/2 * self.r12
        z = 0
        L4 = np.array([x, y, z])
        L5 = np.array([x, -y, z])

        #### L1, L2 and L3 ####
        def f(csi, pi2):
            return (1 - pi2) * (csi + pi2)/np.abs(csi + pi2)**3 + pi2 * (csi + pi2 - 1)/np.abs(csi + pi2 -1)**3 - csi
        
        sign = 1 if self.m1 >= self.m2 else -1
        L1_csi0 = (self.pi1 - self.pi2)
        L2_csi0 = (1 + min(self.pi1, self.pi2)) * sign
        L3_csi0 = (1 + min(self.pi1, self.pi2)) * -sign
        L1_csi = fsolve(f, L1_csi0, args=(self.pi2))[0]
        L1 = np.array([L1_csi * self.r12, 0, 0])
        L2_csi = fsolve(f, L2_csi0, args=(self.pi2))[0]
        L2 = np.array([L2_csi * self.r12, 0, 0])
        L3_csi = fsolve(f, L3_csi0, args=(self.pi2))[0]
        L3 = np.array([L3_csi * self.r12, 0, 0])

        #print(L1_csi0, L1_csi)
        #print(L2_csi0, L2_csi)
        #print(L3_csi0, L3_csi)

        if add_points:
            self.add_point_in_orbital_plane(*self.convert_cartesian_to_polar(L1), name="L1")
            self.add_point_in_orbital_plane(*self.convert_cartesian_to_polar(L2), name="L2")
            self.add_point_in_orbital_plane(*self.convert_cartesian_to_polar(L3), name="L3")
            self.add_point_in_orbital_plane(*self.convert_cartesian_to_polar(L4), name="L4")
            self.add_point_in_orbital_plane(*self.convert_cartesian_to_polar(L5), name="L5")
        return L1, L2, L3, L4, L5

    def jacobi_constant(self, r_vec, v):
        """
        Calculates the Jacobi constant.
        """
        r1 = np.linalg.norm(np.array([r_vec[0]+self.pi2*self.r12, r_vec[1], r_vec[2]]))
        r2 = np.linalg.norm(np.array([r_vec[0]-self.pi1*self.r12, r_vec[1], r_vec[2]]))
        C = (v**2 - self.omega**2 * (r_vec[0]**2 + r_vec[1]**2) - 2*self.mu1/r1 - 2*self.mu2/r2)/2
        return C
    
    def v_for_C(self, r_vec, C):
        """
        Calculates the velocity vector from some specific Jacobi constant (energy).
        """
        r1 = np.linalg.norm(np.array([r_vec[0]+self.pi2*self.r12, r_vec[1], r_vec[2]]))
        r2 = np.linalg.norm(np.array([r_vec[0]-self.pi1*self.r12, r_vec[1], r_vec[2]]))
        v = np.sqrt(2*C + self.omega**2 * (r_vec[0]**2 + r_vec[1]**2) + 2*self.mu1/r1 + 2*self.mu2/r2)
        return v

    def trajectory(self, r_vec_0, v_vec_0, t_span, t_eval=None, method='RK45', add_trajectory_points=True):
        """
        Simulates the trajectory of the spacecraft in the rotating frame.
        """
        def f(t, y):
            y1 = y[0] # x
            y2 = y[1] # y
            y3 = y[2] # z

            y4 = y[3] # x'
            y5 = y[4] # y'
            y6 = y[5] # z'
            
            r1 = np.sqrt((y1 + self.pi2*self.r12)**2 + y2**2 + y3**2)
            r2 = np.sqrt((y1 - self.pi1*self.r12)**2 + y2**2 + y3**2)

            y1_dot = y4
            y2_dot = y5
            y3_dot = y6
            y4_dot = 2*self.omega*y5 + self.omega**2*y1 - self.mu1/r1**3 * (y1 + self.pi2*self.r12) - self.mu2/r2**3 * (y1 - self.pi1*self.r12)
            y5_dot = -2*self.omega*y4 + self.omega**2*y2 - self.mu1/r1**3 * y2 - self.mu2/r2**3 * y2
            y6_dot = - self.mu1/r1**3 * y3 - self.mu2/r2**3 * y3
            return np.concatenate(((y1_dot, y2_dot, y3_dot), (y4_dot, y5_dot, y6_dot)))
        
        sol = solve_ivp(f, t_span, np.concatenate((r_vec_0, v_vec_0)), t_eval=t_eval, method=method)
        if add_trajectory_points:
            # Converter cada ponto (x,y) para (r,theta)
            for point in sol.y[:2].T:  # Pegamos apenas x e y, ignorando z
                r, theta = self.convert_cartesian_to_polar(point)
                self.add_trajectory_points(r, theta)
        return sol
    