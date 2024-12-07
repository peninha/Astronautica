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

        h_vec = np.cross(r_vec, v_vec)
        h = np.linalg.norm(h_vec)
        
        # Vector pointing to the periapsis (eccentricity vector)
        e_vec = np.cross(v_vec, h_vec)/mu - r_vec/np.linalg.norm(r_vec)
        e = np.linalg.norm(e_vec)
        rp = cls.rp_from_mu_h_e(mu, h, e)
        ra = cls.ra_from_mu_h_e(mu, h, e)
        a = cls.a_from_rp_ra(rp, ra)
        epsilon = cls.epsilon_from_mu_a(mu, a)
        p = cls.p_from_mu_h(mu, h)
        b = cls.b_from_a_e(a, e)
        T = cls.T_from_mu_a(mu, a)

        cos_theta = np.dot(e_vec, r_vec)/(e * np.linalg.norm(r_vec))
        sin_theta = np.dot(h_vec, np.cross(e_vec, r_vec))/(e * h * np.linalg.norm(r_vec))
        theta = np.degrees(np.arctan2(sin_theta, cos_theta))

        r = h**2 / (mu) * 1/(1 + e * np.cos(np.radians(theta)))

        orb_params = {
            "mu": mu,
            "m1": m1,
            "m2": m2,
            "h": h,
            "e": e,
            "rp": rp,
            "ra": ra,
            "a": a,
            "b": b,
            "p": p,
            "T": T,
            "epsilon": epsilon,
            "h_vec": h_vec,
            "e_vec": e_vec,
            "r": r,
            "theta": theta
        }
        return orb_params

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

    #####   Geocentric Equatorial Frame   #####
    @staticmethod
    def ra_dec_from_r_vec(r_vec):
        """
        Calculates the right ascension and declination from the position vector.
        """
        r = np.linalg.norm(r_vec)
        ra = np.mod(np.degrees(np.arctan2(r_vec[1], r_vec[0])), 360)
        dec = np.mod(np.degrees(np.arcsin(r_vec[2]/r)), 360)
        return ra, dec

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
            return np.degrees(2 * np.arctanh(np.sqrt((self.e - 1) / (self.e + 1)) * np.tan(theta/2)))
        else:
            raise ValueError("Can't calculate eccentric anomaly from true anomaly for parabolic orbits.")

    def Q_at_delta_t(self, delta_t):
        """
        Finds the universal anomaly at a given time.
        """
        Q0 = np.sqrt(self.mu) * np.abs(self.alpha) * delta_t
        
        def f(Q):
            z = self.alpha * Q**2
            Q = (self.r0*self.vr0 / np.sqrt(self.mu) * Q**2 * self.C(z) + 
                    (1 - self.alpha*self.r0) * Q**3 * self.S(z) + 
                    self.r0 * Q - np.sqrt(self.mu) * delta_t)
            return Q
        def f_dot(Q):
            z = self.alpha * Q**2
            Q_dot = (self.r0*self.vr0 / np.sqrt(self.mu) * Q * (1 - self.alpha *Q**2*self.S(z)) +
                    (1 - self.alpha*self.r0) * Q**2 * self.C(z) + self.r0)
            return Q_dot
        
        Q = root_scalar(f, fprime=f_dot, x0=Q0, method='newton').root
        return Q


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

    def lagrange_coefficients_universal(self, r_vec_0, v_vec_0, Q, delta_t):
        """
        Calculates the Lagrange coefficients using universal anomaly.

        :param r_vec_0: Initial position vector (km)
        :param v_vec_0: Initial velocity vector (km/s)
        :param Q: Universal anomaly (km^1/2)
        :param delta_t: Time variation (s)
        :return: f, g, f_dot, g_dot
        """
        r0 = np.linalg.norm(r_vec_0)
        v0 = np.linalg.norm(v_vec_0)
        vr0 = np.dot(r_vec_0, v_vec_0) / r0
        alpha = self.alpha_from_mu_r_v(self.mu, r0, v0)
        z = alpha * Q**2

        f = 1 - Q**2/r0 * self.C(z)
        g = delta_t - Q**3/np.sqrt(self.mu) * self.S(z)
        
        r_vec = f*r_vec_0 + g*v_vec_0
        r = np.linalg.norm(r_vec)

        f_dot = -np.sqrt(self.mu)/(r*r0) * (alpha*Q**3 * self.S(z) - Q)
        g_dot = 1 - Q**2/r * self.C(z)

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
            M = self.M_from_t(t)
            E = self.E_from_M(M)
            return np.mod(self.theta_from_E(E), 360)
        else: # parabolic
            M = self.M_from_t(t)
            M = np.radians(M)
            z = (3*M + np.sqrt(1 + 9*M**2))**(1/3)
            theta = np.degrees(2 * np.arctan(z - 1/z))
            return np.mod(theta, 360)


    #######   Time evolution   #######
    def t_at_M(self, M):
        """
        Calculates the time using the mean anomaly.
        """
        M = np.radians(M)
        if self.e < 1:
            return M * self.T / (2 * np.pi)
        elif self.e == 1:
            return M * self.h**3 / self.mu**2
        else:
            return M * self.h**3 / (self.mu**2 * (self.e**2 - 1)**(3/2))
    
    def t_at_theta(self, theta):
        """
        Calculates the time using the true anomaly.
        """
        if self.e < 1 or self.e > 1: # elliptical or hyperbolic
            E = self.E_from_theta(theta)
            M = self.M_from_E(E)
            return self.t_from_M(M)
        else: # parabolic
            M = self.M_from_theta(theta)
            return self.t_from_M(M)


    #######   State vectors   #######
    def r_vec_at_theta(self, theta, frame="perifocal"):
        """
        Calculates the position vector at a given true anomaly, in a given frame of reference.
        """
        theta = np.radians(theta)
        if frame == "perifocal":
            return self.h**2 / (self.mu * (1 + self.e*np.cos(theta))) * np.array([np.cos(theta), np.sin(theta), 0])
        elif frame == "geocentric equatorial":
            raise NotImplementedError("Not implemented yet.")
        else:
            raise ValueError("Invalid frame of reference.")
    
    def v_vec_at_theta(self, theta, frame="perifocal"):
        """
        Calculates the velocity vector at a given true anomaly, in a given frame of reference.
        """
        theta = np.radians(theta)
        if frame == "perifocal":
            return self.mu / self.h * np.array([-np.sin(theta), self.e + np.cos(theta), 0])
        elif frame == "geocentric equatorial":
            raise NotImplementedError("Not implemented yet.")
        else:
            raise ValueError("Invalid frame of reference.")    

    def state_at_Q(self, Q, delta_t):
        """
        Finds the position and velocity vectors at a given universal anomaly and delta time.
        """
        z = self.alpha * Q**2
        
        ## Lagrange coefficients ##
        f = 1 - Q**2/self.r0 * self.C(z)
        g = delta_t - Q**3/np.sqrt(self.mu) * self.S(z)
        
        r_vec = f*self.r_vec_0 + g*self.v_vec_0
        r = np.linalg.norm(r_vec)

        ## Lagrange coefficients derivatives ##
        f_dot = -np.sqrt(self.mu)/(r*self.r0) * (self.alpha*Q**3 * self.S(z) - Q)
        g_dot = 1 - Q**2/r * self.C(z)

        v_vec = f_dot*self.r_vec_0 + g_dot*self.v_vec_0

        return r_vec, v_vec

    def state_at_theta(self, theta, frame="perifocal"):
        """
        Calculates the state vector at a given true anomaly, in a given frame of reference.
        """
        if frame == "perifocal":
            r = self.r_vec_at_theta(theta, frame)
            v = self.v_vec_at_theta(theta, frame)
        elif frame == "geocentric equatorial":
            raise NotImplementedError("Not implemented yet.")
        else:
            raise ValueError("Invalid frame of reference.")
        return r, v

    def state_at_t(self, t):
        """
        Finds the position and velocity vectors at a given time.
        """
        delta_t = t - self.t0
        Q = self.Q_at_delta_t(delta_t)
        r_vec, v_vec = self.state_at_Q(Q, delta_t)
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
        Q = self.Q_at_delta_t(delta_t)
        return self.state_at_Q(Q, delta_t)


    #######   Plotting   #######
    def plot(self, plot_points=True, plot_positions=True, plot_velocities=False):
        """
        Plots the orbit with optional points and velocities.

        :param plot_points: Boolean to plot the points
        :param plot_positions: Boolean to plot the already added orbital positions
        :param plot_velocities: Boolean to plot the velocities
        """
        # Criar array de anomalia verdadeira
        if self.e == 1:
            theta_min = -120
            theta_max = 120
            if plot_positions:
                for pos in self.positions:
                    if pos['theta'] < theta_min:
                        theta_min = pos['theta']
                    if pos['theta'] > theta_max:
                        theta_max = pos['theta']
            theta = np.linspace(theta_min, theta_max, 1000)
        elif self.e > 1:
            epsilon = 15
            theta_min = -self.theta_inf() + epsilon
            theta_max = self.theta_inf() - epsilon
            if plot_positions:
                for pos in self.positions:
                    if pos['theta'] < theta_min:
                        theta_min = pos['theta']
                    if pos['theta'] > theta_max:
                        theta_max = pos['theta']
            theta = np.linspace(theta_min, theta_max, 1000)
        else:
            theta = np.linspace(0, 360, 1000)
        
        # Calcular raio para cada ângulo
        r = self.r_at_theta(theta)
        
        # Converter para coordenadas cartesianas
        x = r * np.cos(np.radians(theta))
        y = r * np.sin(np.radians(theta))
        
        # Criar o plot
        plt.figure(figsize=(8, 8))
        plt.plot(x, y, 'b-', label='Orbit')
        
        # Plotar o corpo central com o raio correto se fornecido
        if self.body1radius is not None:
            circle = plt.Circle((0, 0), self.body1radius, color='r', alpha=0.3, label='Central Body')
            plt.gca().add_patch(circle)
        else:
            plt.plot(0, 0, 'ro', label='Corpo Central')
    
        # list of colors
        cores = ['red', 'orange', 'purple', 'brown', 'deeppink', 'gray', 'olive', 'deepskyblue', 'blue', 'green']
                
        if plot_points:
            for i, point in enumerate(self.points_in_orbital_plane):
                r = point['r']
                theta = np.radians(point['theta'])
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                label = point.get('name', f'Point {i+1}')
                cor = 'black'
                plt.plot(x, y, 'o', color=cor, label=label)

        if plot_positions:
            for i, position in enumerate(self.positions):
                r = position['r']
                theta = np.radians(position['theta'])
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                label = position.get('name', f'Position {i+1}')
                if position['t'] is not None:
                    label += f' (t={position["t"]:.2f}s)'
                cor = cores[i % len(cores)]  # Recycle colors if there are more points than colors
                # Plot the secondary body if provided
                if self.body2radius is not None:
                    circle2 = plt.Circle((x, y), self.body2radius, color=cor, alpha=0.3, label=label)
                    plt.gca().add_patch(circle2)
                else:
                    plt.plot(x, y, 'o', color=cor, label=label)
                
                # Plot velocity vector
                if plot_velocities:
                    v = position['v']
                    v_vec = self.v_vec_at_theta(position['theta'])
                    vx = v_vec[0]
                    vy = v_vec[1]
                    # Normalize the vector for better visualization
                    scale = self.p/20  # Velocity vector scale
                    plt.arrow(x, y, vx*scale, vy*scale, 
                            color=cor, width=0.1, head_width=1*scale, 
                            head_length=1.5*scale, alpha=0.8)

        plt.grid(True)
        plt.axis('equal')
        plt.xlabel('x (km)')
        plt.ylabel('y (km)')
        plt.title('Orbit')
        plt.legend()


    ####### String representation #######
    def __str__(self):
        """
        Returns a string representation of the orbital parameters.
        """
        # Criar uma lista com todos os atributos da classe
        atributos = vars(self)
        
        # String inicial
        resultado = "Parâmetros Orbitais:\n"
        
        # Adicionar cada atributo à string
        for nome, valor in atributos.items():
            # Pular atributos que são None
            if valor is None:
                continue
                
            # Formatar números como float com 3 casas decimais
            if isinstance(valor, (int, float)):
                resultado += f"{nome}: {valor:.3f}\n"
            else:
                resultado += f"{nome}: {valor}\n"
                
        return resultado
    



class Orbit(OrbitalBase):
    def __init__(self, mu=None, m1=None, m2=0, a=None, e=None, rp=None, ra=None, h=None, body1radius=None, body2radius=None):
        """
        Initializes the Orbit class with primary orbital parameters.
        :param mu: Standard gravitational parameter (km^3/s^2)
        :param a: Semi-major axis (km)
        :param e: Eccentricity (dimensionless)
        :param rp: Periapsis distance (km)
        :param ra: Apoapsis distance (km)
        :param h: Specific angular momentum (km^2/s)
        :param epsilon: Specific orbital energy (km^2/s^2)
        :param body1radius: Radius of the central body (km)
        :param body2radius: Radius of the orbiting body (km)
        """
        super().__init__()
        self.body1radius = body1radius
        self.body2radius = body2radius

        if mu is not None:
            self.mu = mu  # Gravitational parameter
        elif m1 is not None:
            self.mu = self.mu_from_m1_m2(m1, m2)
        else:
            raise ValueError("Provide either mu or m1 (and optionally m2) to calculate mu.")
        self.m1 = m1
        self.m2 = m2
        
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
            self.e = self.e_from_mu_h_rp(self.mu, self.h, self.rp)
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
    def init_from_state_vectors(cls, r_vec, v_vec, mu=None, m1=None, m2=0, name=None, body1radius=None, body2radius=None):
        """
        Given a position and velocity vectors represented in some inertial frame of reference
        centered in the primary body, creates a new instance of Orbit.

        :param r_vec: Position vector in the inertial frame of reference (km)
        :param v_vec: Velocity vector in the inertial frame of reference (km/s)
        :param mu: Standard gravitational parameter (km³/s²)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param name: Name of the body
        :param body1radius: Radius of the central body (km) (optional)
        :param body2radius: Radius of the orbiting body (km) (optional)
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Provide mu or m1 (and optionally m2) to calculate mu.")
        
        orbital_parameters = cls.orbital_parameters_from_state_vectors(r_vec, v_vec, mu)
        e = orbital_parameters['e']
        h = orbital_parameters['h']
        theta = orbital_parameters['theta']

        orbit_instance = cls(mu=mu, m1=m1, m2=m2, e=e, h=h, body1radius=body1radius, body2radius=body2radius)
        orbit_instance.add_orbital_position(theta=theta, name="Position 1")
        
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

        e = Orbit.e_from_mu_r_v_theta(mu, r, v, theta)
        h = Orbit.h_from_mu_r_e_theta(mu, r, e, theta)

        orbit_instance = cls(mu=mu, m1=m1, m2=m2, e=e, h=h, body1radius=body1radius, body2radius=body2radius)
        orbit_instance.add_orbital_position(theta=theta, name="Position 1")
        
        return orbit_instance

    ############# Adding Points #############
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
        self.body1radius = body1radius
        self.body2radius = body2radius
        self.omega = np.sqrt(self.mu/self.r12**3)

    def lagrange_points(self, plot=False):
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

        print(L1_csi0, L1_csi)
        print(L2_csi0, L2_csi)
        print(L3_csi0, L3_csi)

        if plot:
            self.plot(points=[L1, L2, L3, L4, L5])
        
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

    def trajectory(self, r_vec_0, v_vec_0, t_span, t_eval=None, method='RK45', plot=False):
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
        if plot:
            self.plot(points=sol.y[:3].T)
        return sol
    