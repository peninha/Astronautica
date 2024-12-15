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
    def orbital_parameters_from_state_vectors(cls, r_vec_bc, v_vec_bc, mu=None, m1=None, m2=0):
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
        h_vec_bc = cls.h_vec_from_state_vectors(r_vec_bc, v_vec_bc)
        h = np.linalg.norm(h_vec_bc)

        # Inclination   
        i = cls.i_from_h_vec(h_vec_bc)

        # Node vector
        n_vec_bc = cls.n_vec_from_h_vec(h_vec_bc)
        n = np.linalg.norm(n_vec_bc)

        # Right ascension of the ascending node
        Omega = cls.Omega_from_n_vec(n_vec_bc)

        # Eccentricity vector
        e_vec_bc = cls.e_vec_from_mu_state_vectors(mu, r_vec_bc, v_vec_bc)
        e = np.linalg.norm(e_vec_bc)

        # Argument of periapsis
        omega = cls.omega_from_e_vec_n_vec(e_vec_bc, n_vec_bc)

        # True anomaly
        theta0 = cls.theta_from_h_vec_e_vec_r_vec(h_vec_bc, e_vec_bc, r_vec_bc)

        # Distance
        r = np.linalg.norm(r_vec_bc)

        # Velocity
        v = np.linalg.norm(v_vec_bc)

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
            "theta0": theta0,
            "h_vec_bc": h_vec_bc,
            "e_vec_bc": e_vec_bc,
            "n_vec_bc": n_vec_bc,
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
    def i_from_h_vec(h_vec_bc):
        """
        Calculates the inclination using specific angular momentum vector.
        h vector must be in the bodycentric frame.
        """
        return np.degrees(np.arccos(h_vec_bc[2]/np.linalg.norm(h_vec_bc)))

    #####   N   #####
    @staticmethod
    def n_vec_from_h_vec(h_vec_bc):
        """
        Calculates the node vector using specific angular momentum vector.
        h vector must be in the bodycentric frame.
        """
        return np.cross(np.array([0, 0, 1]), h_vec_bc)

    #####   Omega   #####
    @staticmethod
    def Omega_from_n_vec(n_vec_bc):
        """
        Calculates the right ascension of the ascending node using node vector.
        Node vector must be in the bodycentric frame.
        """
        n = np.linalg.norm(n_vec_bc)
        if n == 0:
            Omega = 0
        else:
            cos_Omega = n_vec_bc[0]/n
            sin_Omega = n_vec_bc[1]/n
            Omega = np.mod(np.degrees(np.arctan2(sin_Omega, cos_Omega)), 360)
        return Omega

    #####   omega   #####
    @staticmethod
    def omega_from_e_vec_n_vec(e_vec, n_vec):
        """
        Calculates the argument of periapsis using eccentricity vector and node vector.
        Both vectors must be in the same reference frame.
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
    def theta_from_h_vec_e_vec_r_vec(h_vec, e_vec, r_vec):
        """
        Calculates the true anomaly using specific angular momentum vector, eccentricity vector and position vector.
        All vectors must be in the same reference frame.
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
        Both vectors must be in the same reference frame.
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
    def h_from_state_vectors(r_vec, v_vec):
        """
        Calculates the specific angular momentum using position and velocity vectors.
        """
        return np.linalg.norm(OrbitalBase.h_vec_from_state_vectors(r_vec, v_vec))
    
    @staticmethod
    def h_vec_from_state_vectors(r_vec, v_vec):
        """
        Calculates the specific angular momentum vector using position and velocity vectors.
        Both vectors must be in the same reference frame.
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

    @staticmethod
    def omega_dot_oblateness_correction(J2, mu, R_body, i, e, a):
        """
        Calculates the oblateness correction factor using gravitational parameter, distance, velocity and body rotation rate.
        """
        if e >= 1:
            raise ValueError("Can't calculate oblateness correction for parabolic or hyperbolic orbits.")
        i = np.radians(i)
        coeficient = - 3/2 * np.sqrt(mu) * J2 * R_body**2 / (a**(7/2) * (1 - e**2)**2)
        Omega_dot = np.degrees(coeficient * np.cos(i))
        omega_dot = np.degrees(coeficient * (5/2 * np.sin(i)**2 - 2)) 
        return Omega_dot, omega_dot

    @staticmethod
    def step_function(var, var_period):
        """
        Step function for variable correction of periodic variables outside
        the range of -var_period/2 and var_period/2.
        Returns:
        -var_period/2 <= var <= var_period/2: 0
        -3*var_period/2 <= var < -var_period/2: -1
        -5*var_period/2 <= var < -3*var_period/2: -2
        ...
        var_period/2 < var <= 3*var_period/2: 1
        5*var_period/2 < var <= 9*var_period/2: 2
        ...
        """
        if -var_period/2 <= var <= var_period/2:
            return 0
        if var < var_period/2:
            return (var + var_period/2) // var_period
        else:
            return np.abs((-var + var_period/2) // var_period)

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

    @classmethod
    #######   Lagrange coefficients   #######
    def lagrange_coefficients_true(cls, mu, r_vec_0, v_vec_0, delta_theta):
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

        r = h**2 / (mu * (1 + (h**2/(mu*r0) - 1)*np.cos(delta_theta) - h*vr0/mu * np.sin(delta_theta)))
        f = 1 - mu*r/h**2 * (1 - np.cos(delta_theta))
        g = r*r0/h * np.sin(delta_theta)
        f_dot = mu/h * (1 - np.cos(delta_theta))/np.sin(delta_theta) * (mu/h**2 * (1 - np.cos(delta_theta)) - 1/r0 - 1/r)
        g_dot = 1 - mu*r0/h**2 * (1 - np.cos(delta_theta))

        return f, g, f_dot, g_dot

    @classmethod
    def lagrange_coefficients_universal(cls, mu, r_vec_0, v_vec_0, Qui, delta_t):
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
        alpha = cls.alpha_from_mu_r_v(mu, r0, v0)
        z = alpha * Qui**2

        f = 1 - Qui**2/r0 * cls.C(z)
        g = delta_t - Qui**3/np.sqrt(mu) * cls.S(z)
        
        r_vec = f*r_vec_0 + g*v_vec_0
        r = np.linalg.norm(r_vec)

        f_dot = np.sqrt(mu)/(r*r0) * (alpha*Qui**3 * cls.S(z) - Qui)
        g_dot = 1 - Qui**2/r * cls.C(z)

        return f, g, f_dot, g_dot

    @classmethod
    def Qui_at_delta_t(cls, r0, vr0, alpha, delta_t, mu=None, m1=None, m2=0):
        """
        Finds the universal anomaly at a given time.
        """
        if mu is None:
            mu = cls.mu_from_m1_m2(m1, m2)
        elif m1 is None or m2 is None:
            raise ValueError("Can't calculate universal anomaly without gravitational parameter or masses.")
        
        Q0 = np.sqrt(mu) * np.abs(alpha) * delta_t
        
        def f(Qui):
            z = alpha * Qui**2
            Qui = (r0*vr0 / np.sqrt(mu) * Qui**2 * cls.C(z) + 
                    (1 - alpha*r0) * Qui**3 * cls.S(z) + 
                    r0 * Qui - np.sqrt(mu) * delta_t)
            return Qui
        def f_dot(Qui):
            z = alpha * Qui**2
            Qui_dot = (r0*vr0 / np.sqrt(mu) * Qui * (1 - alpha * Qui**2 * cls.S(z)) +
                    (1 - alpha*r0) * Qui**2 * cls.C(z) + r0)
            return Qui_dot
        
        Qui = root_scalar(f, fprime=f_dot, x0=Q0, method='newton').root
        return Qui
    
    @classmethod
    def state_vectors_at_Qui(cls, r_vec_0, v_vec_0, Qui, delta_t, mu=None, m1=None, m2=0):
        """
        Finds the position and velocity vectors at a given universal anomaly and delta time.
        State vectors are in the same reference frame as the initial state vectors, which is perifocal.

        :param r_vec_0: Initial position vector (km)
        :param v_vec_0: Initial velocity vector (km/s)
        :param alpha: Alpha parameter, inverse of the semi-major axis (1/km)
        :param Qui: Universal anomaly (km^1/2)
        :param delta_t: Time change (s)
        :return: Position and velocity vectors at the given universal anomaly and delta time (km and km/s)
        """
        if mu is None:
            mu = cls.mu_from_m1_m2(m1, m2)
        elif m1 is None or m2 is None:
            raise ValueError("Can't calculate state vectors at universal anomaly without gravitational parameter or masses.")
        
        alpha = cls.alpha_from_mu_r_v(mu, r_vec_0, v_vec_0)
        z = alpha * Qui**2
        r0 = np.linalg.norm(r_vec_0)

        ## Lagrange coefficients ##
        f = 1 - Qui**2/r0 * cls.C(z)
        g = delta_t - Qui**3/np.sqrt(mu) * cls.S(z)
        
        r_vec = f*r_vec_0 + g*v_vec_0
        r = np.linalg.norm(r_vec)

        ## Lagrange coefficients derivatives ##
        f_dot = np.sqrt(mu)/(r*r0) * (alpha*Qui**3 * cls.S(z) - Qui)
        g_dot = 1 - Qui**2/r * cls.C(z)

        v_vec = f_dot*r_vec_0 + g_dot*v_vec_0

        return r_vec, v_vec
    
    @classmethod
    def state_vectores_after_delta_theta(cls, r_vec_0, v_vec_0, delta_theta, mu=None, m1=None, m2=0):
        """
        Calculates the position and velocity vectors after a given true anomaly change.
        Useful for propagating orbits in the original frame of reference, regardless of the reference frame used.

        :param r_vec_0: Initial position vector (km)
        :param v_vec_0: Initial velocity vector (km/s)
        :param delta_theta: Change in true anomaly (degrees)
        :return: Position and velocity vectors after the change (km and km/s)
        """
        if mu is None:
            mu = cls.mu_from_m1_m2(m1, m2)
        elif m1 is None or m2 is None:
            raise ValueError("Can't calculate state vectors after true anomaly change without gravitational parameter or masses.")
        
        f, g, f_dot, g_dot = cls.lagrange_coefficients_true(mu, r_vec_0, v_vec_0, delta_theta)
        r_vec = f*r_vec_0 + g*v_vec_0
        v_vec = f_dot*r_vec_0 + g_dot*v_vec_0
        return r_vec, v_vec
    
    @classmethod
    def state_vectores_after_delta_t(cls, r_vec_0, v_vec_0, delta_t, mu=None, m1=None, m2=0):
        """
        Calculates the position and velocity vectors after a given time change.
        Useful for propagating orbits in the original frame of reference, regardless of the reference frame used.

        :param r_vec_0: Initial position vector (km)
        :param v_vec_0: Initial velocity vector (km/s)
        :param delta_t: Time change (s)
        :return: Position and velocity vectors after the change (km and km/s)
        """
        if mu is None:
            mu = cls.mu_from_m1_m2(m1, m2)
        elif m1 is None or m2 is None:
            raise ValueError("Can't calculate state vectors after time change without gravitational parameter or masses.")
        
        alpha = cls.alpha_from_mu_r_v(mu, r_vec_0, v_vec_0)
        r0 = np.linalg.norm(r_vec_0)
        vr0 = np.dot(r_vec_0, v_vec_0) / r0
        Qui = cls.Qui_at_delta_t(r0, vr0, alpha, delta_t, mu=mu)
        return cls.state_vectors_at_Qui(r_vec_0, v_vec_0, Qui, delta_t, mu=mu)

    #####   Conversions   #####
    @staticmethod
    def color_gradient(value, colors=[(1,0.5,0), (0,0,1)]):
        """
        Returns an interpolated color based on a number between 0 and 1.

        :param colors: List of RGB colors to interpolate (each color as tuple of 3 values 0-1) 
        :param value: Value between 0 and 1 to interpolate the color
        :return: Interpolated RGB color tuple
        """
        if value < 0 or value > 1:
            raise ValueError("Value must be between 0 and 1.")
            
        if len(colors) < 2:
            raise ValueError("At least 2 colors are needed for interpolation.")
            
        # Find the two colors to interpolate between
        num_segments = len(colors) - 1
        segment = value * num_segments
        i = int(segment)
        
        # Handle last segment edge case
        if i >= num_segments:
            return colors[-1]
            
        # Interpolate between the two colors
        t = segment - i
        c1 = np.array(colors[i])
        c2 = np.array(colors[i+1])
        
        return tuple(c1 * (1-t) + c2 * t)
    
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
        ra = np.degrees(np.arctan2(r_vec_xyz[1], r_vec_xyz[0]))
        dec = np.degrees(np.arcsin(r_vec_xyz[2]/r))
        return ra, dec

    @staticmethod
    def convert_bc_to_rbc(vec_bc, omega_body, delta_t):
        """
        Converts a body-centric vector to a rotating body-centric vector.

        :param vec_bc: Body-centric vector (km)
        :param omega_body: Body rotation rate (deg/s)
        :param delta_t: Time change (s)
        :return: Rotating body-centric vector (km)
        """
        Q_bc_rbc = OrbitalBase.Q_from_Euler_angles(omega_body * delta_t, pattern="z")
        return Q_bc_rbc @ vec_bc

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
        Calculates the flight path angle from state vectors.
        
        :param r_vec: Position vector (km)
        :param v_vec: Velocity vector (km/s)
        :return: Flight path angle (degrees)
        """
        # Calculate the magnitudes of the vectors
        r = np.linalg.norm(r_vec)
        v = np.linalg.norm(v_vec)
        
        # Calculate the dot product
        r_dot_v = np.dot(r_vec, v_vec)
        
        # Calculate the angle using arccos of the normalized dot product
        gamma = np.arccos(r_dot_v/(r*v))
        
        # Convert to degrees and subtract from 90° to obtain the flight path angle
        return 90 - np.degrees(gamma)

    def add_oblateness_correction(self, J2=1.0826e-3, body1radius=None):
        """
        Adds the oblateness correction parameters to the orbit.
        """
        if body1radius is None and self.body1radius is None:
            raise ValueError("Can't add oblateness correction without body radius.")
        elif body1radius is not None:
            self.body1radius = body1radius
        
        self.Omega_dot, self.omega_dot = self.omega_dot_oblateness_correction(J2, self.mu, self.body1radius, self.i, self.e, self.a)
        self.J2 = J2

    def add_main_body_rotation(self, omega_body=np.degrees(7.292115e-5)):
        """
        Adds the main body rotation to the orbit.
        """
        self.omega_body = omega_body

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

    def M_at_t_orbit(self, t_orbit):
        """
        Calculates the mean anomaly using the orbit time.
        t_orbit=0 is the time of the perigee passage.
        """
        if self.e < 1:
            return np.degrees(2 * np.pi * t_orbit / self.T)
        elif self.e == 1:
            return np.degrees(self.mu**2 * t_orbit/ self.h**3)
        else:
            return np.degrees(self.mu**2 / self.h**3 * (self.e**2 - 1)**(3/2) * t_orbit)

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

    #######   Rotation   #######
    def Q_bc_peri(self, t_clock):
        """
        Calculates the direction cosine matrix from perifocal frame to bodycentric frame at time t0.
        """
        delta_t = t_clock - self.t0_clock
        return self.Q_from_Euler_angles(self.Omega + self.Omega_dot * delta_t,
                                        self.i, 
                                        self.omega + self.omega_dot * delta_t,
                                        pattern="classical")

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

    def theta_at_t_clock(self, t_clock):
        """
        Calculates the true anomaly at a given clock time.
        """
        delta_t = t_clock - self.t0_clock
        t_orbit = self.t0_orbit + delta_t
        if self.e < 1 or self.e > 1: # elliptical or hyperbolic
            M = self.M_at_t_orbit(t_orbit)
            E = self.E_from_M(M)
            return self.theta_from_E(E) + self.step_function(t_orbit, self.T) * 360
        else: # parabolic
            M = self.M_at_t_orbit(t_orbit)
            M = np.radians(M)
            z = (3*M + np.sqrt(1 + 9*M**2))**(1/3)
            return np.degrees(2 * np.arctan(z - 1/z))

    def theta_at_state_vectors(self, r_vec, v_vec):
        """
        Calculates the true anomaly angle from state vectors,
        regardless of the reference frame used.
        """
        h_vec = self.h_vec_from_state_vectors(r_vec, v_vec)
        h = np.linalg.norm(h_vec)
        e_vec = self.e_vec_from_mu_state_vectors(self.mu, r_vec, v_vec)
        e = np.linalg.norm(e_vec)
        if e == 0:
            theta = 0
        else:
            cos_theta = np.dot(e_vec, r_vec)/(e * np.linalg.norm(r_vec))
            sin_theta = np.dot(h_vec, np.cross(e_vec, r_vec))/(e * h * np.linalg.norm(r_vec))
            theta = np.mod(np.degrees(np.arctan2(sin_theta, cos_theta)), 360)
        return theta


    #######   Time evolution   #######
    def t_orbit_at_M(self, M):
        """
        Calculates the orbit time using the mean anomaly.
        t_orbit=0 is the time of the perigee passage.
        """
        M = np.radians(M)
        if self.e < 1:
            t = M * self.T / (2 * np.pi)
        elif self.e == 1:
            t = M * self.h**3 / self.mu**2
        else:
            t = M * self.h**3 / (self.mu**2 * (self.e**2 - 1)**(3/2))
        return t
    
    def t_orbit_at_theta(self, theta):
        """
        Calculates the orbit time using the true anomaly.
        t_orbit=0 is the time of the perigee passage.
        """
        if self.e < 1 or self.e > 1: # elliptical or hyperbolic
            E = self.E_from_theta(theta)
            M = self.M_from_E(E)
            return self.t_orbit_at_M(M) + self.step_function(theta, 360) * self.T
        else: # parabolic
            M = self.M_from_theta(theta)
            return self.t_orbit_at_M(M)

    def t_clock_at_theta(self, theta):
        """
        Calculates the clock time using the true anomaly.
        """
        t_orbit = self.t_orbit_at_theta(theta)
        return t_orbit + self.time_offset

    #######   State vectors   #######
    def r_vec_at_theta(self, theta, frame="perifocal"):
        """
        Calculates the position vector at a given true anomaly, in a given frame of reference.
        """
        r_vec, v_vec = self.state_vectors_at_theta(theta, frame)
        return r_vec
    
    def v_vec_at_theta(self, theta, frame="perifocal"):
        """
        Calculates the velocity vector at a given true anomaly, in a given frame of reference.
        """
        r_vec, v_vec = self.state_vectors_at_theta(theta, frame)
        return v_vec

    def state_vectors_at_theta(self, theta, frame="perifocal"):
        """
        Calculates the state vectors at a given true anomaly, in a given frame of reference, at time t0.

        :param theta: True anomaly (degrees)
        :param frame: Frame of reference to calculate the position vector
        """
        if frame == "perifocal":
            theta_rad = np.radians(theta)
            r_vec = self.h**2 / (self.mu * (1 + self.e*np.cos(theta_rad))) * np.array([np.cos(theta_rad), np.sin(theta_rad), 0])
            v_vec = self.mu / self.h * np.array([-np.sin(theta_rad), self.e + np.cos(theta_rad), 0])
        elif frame == "perifocal_t0":
            if not self.omega_dot and not self.Omega_dot:
                return self.state_vectors_at_theta(theta, frame="perifocal")
            else:
                r_vec_bc, v_vec_bc = self.state_vectors_at_theta(theta, frame="bodycentric")
                r_vec = self.Q_bc_perit0 @ r_vec_bc
                v_vec = self.Q_bc_perit0 @ v_vec_bc
        elif frame == "bodycentric":
            if not self.omega_dot and not self.Omega_dot:
                Q_bc_peri = self.Q_bc_perit0
            else:
                t_orbit = self.t_orbit_at_theta(theta)
                t_clock = t_orbit + self.time_offset
                Q_bc_peri = self.Q_bc_peri(t_clock)
            r_vec_peri, v_vec_peri = self.state_vectors_at_theta(theta, frame="perifocal")
            r_vec = Q_bc_peri.T @ r_vec_peri
            v_vec = Q_bc_peri.T @ v_vec_peri

        elif frame == "rotatingBodycentric":
            t_orbit = self.t_orbit_at_theta(theta)
            t_clock = t_orbit + self.time_offset
            delta_t = t_clock - self.t0_clock
            r_vec = self.convert_bc_to_rbc(self.r_vec_at_theta(theta, frame="bodycentric"), self.omega_body, delta_t)
            v_vec = self.convert_bc_to_rbc(self.v_vec_at_theta(theta, frame="bodycentric"), self.omega_body, delta_t)
        else:
            raise ValueError("Invalid frame of reference.")
        return r_vec, v_vec

    def state_vectors_at_t_clock(self, t_clock, frame="perifocal"):
        """
        Finds the position and velocity vectors at a given time.
        """
        theta = self.theta_at_t_clock(t_clock)
        return self.state_vectors_at_theta(theta, frame=frame)
    
    #######   Adding Points   #######
    def add_orbital_position(self, theta, t_clock=None, name=None):
        """
        Adds an orbital position to the orbit at a given true anomaly and (optionally) time.

        :param theta: True anomaly (degrees)
        :param t_clock: Time (s), default is None
        :param name: Name of the position, default is None
        """
        r = self.r_at_theta(theta)
        v = self.v_at_r(r)
        if t_clock is None:
            t_clock = self.t_clock_at_theta(theta)
        position = {
            'name': name,
            'r': r,
            'v': v, 
            'theta': theta,
            't_clock': t_clock
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

    def add_zero_state_from_theta(self, theta0=0, t0_clock=0, add_position=True):
        """
        Adds the initial state vectors for trajectory calculations.
        State vectors are always given in perifocal frame.
        """
        self.theta0 = theta0
        self.t0_clock = t0_clock
        self.t0_orbit = self.t_orbit_at_theta(self.theta0)
        self.time_offset = self.t0_clock - self.t0_orbit
        self.r_vec_0, self.v_vec_0 = self.state_vectors_at_theta(theta0, frame="perifocal")
        self.r0 = np.linalg.norm(self.r_vec_0)
        self.v0 = np.linalg.norm(self.v_vec_0)
        self.vr0 = np.dot(self.r_vec_0, self.v_vec_0)/self.r0
        if add_position:
            self.add_orbital_position(self.theta0, self.t0_clock, name="Initial Position")

    def add_zero_state_from_vectors(self, r_vec, v_vec, t0_clock=0, add_position=True):
        """
        Adds the initial state vectors for trajectory calculations.

        :param r_vec: Initial position vector (km)
        :param v_vec: Initial velocity vector (km/s)
        :param t0_clock: Time (s)
        :param add_position: Boolean to add the initial position to the orbit
        """
        theta0 = self.theta_at_state_vectors(r_vec, v_vec)
        self.add_zero_state_from_theta(theta0, t0_clock, add_position)

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
    
    def trajectory(self, theta1=360, t1_clock=None, n_points=20):
        """
        Calculates the trajectory between two true anomalies or times.

        :param theta_1: Final true anomaly (degrees) (optional, default is 360)
        :param t1: Final time (s) - if provided, overrides theta_1 (optional, default is None)
        :param n_points: Number of points to evaluate the trajectory (optional, default is 20)
        :param add_trajectory_points: Boolean to add the trajectory points to the orbit (optional, default is True)
        """
        # Minimum interval (s) between points
        if t1_clock is None:
            # If t1 is not provided, calculate from theta1
            self.t1_orbit = self.t_orbit_at_theta(theta1)
            self.t1_clock = self.time_offset + self.t1_orbit
        else:
            self.t1_clock = t1_clock
            self.t1_orbit = t1_clock - self.time_offset
        
        # Add final position
        theta1 = self.theta_at_t_clock(self.t1_clock)
        self.add_orbital_position(theta1, self.t1_clock, name='Final Position')
        
        # Ensure interval between points is not smaller than min_interval
        min_interval = 10
        n_points = min(n_points, max(2, int((self.t1_clock - self.t0_clock) / min_interval)))
        t_eval = np.linspace(self.t0_clock, self.t1_clock, n_points)

        # Calculate trajectory for each evaluation point
        for t_clock in t_eval:
            theta = self.theta_at_t_clock(t_clock)
            r = self.r_at_theta(theta)
            self.add_trajectory_points(r, theta, t_clock)

    #######   Plotting   #######
    def plot(self, orbit=True, points=True, positions=True, velocities=True, trajectory=True, frame="perifocal", plot3d=True, groundtrack=False):
        """
        Plots the orbit with optional points and velocities.

        :param orbit: Boolean to plot the orbit
        :param points: Boolean to plot the points
        :param positions: Boolean to plot the already added orbital positions
        :param velocities: Boolean to plot the velocities
        :param trajectory: Boolean to plot the trajectory
        :param frame: Frame of reference to plot the orbit

        Frames of reference:
        Perifocal frame is the frame with the origin at the primary body center, with:
            x-axis is along the periapsis
            y-axis is inside the orbital plane
            z-axis is along the angular momentum vector
            if the orbit is precessing, the x-axis continuously tracks the periapsis position,
            making this a non-inertial reference frame
        
        Perifocal_t0 frame is a inertial frame with the origin at the primary body center, with:
            x-axis is along the periapsis at time t0
            y-axis is inside the orbital plane
            z-axis is along the angular momentum vector
            This is an inertial frame. If the orbit is precessing,
            the x-axis will not track the periapsis position, but will remain fixed at the initial periapsis position.

        Bodycentric frame is a inertial frame with the origin at the primary body center, with:
            x-axis is inside the equatorial plane, along a predifined direction. If main body is Earth, vernal equinox is used.
            y-axis is perpendicular to the x-axis and in the equatorial plane
            z-axis is pointing to the north pole of the main body
        
        rotatingBodycentric frame is the same as bodycentric frame, but the x-axis is rotating with the main body, with:
            x-axis is inside the equatorial plane, along a predifined direction. At time t0, the x-axis coincides with bodycentric x-axis.
            y-axis is perpendicular to the x-axis and in the equatorial plane
            z-axis is pointing to the north pole of the main body
            The x-axis will rotate with the main body, making this a non-inertial reference frame.

        RotatingBarycentric frame is non-inertial frame with the origin at the barycenter of the system, with:
            x-axis is pointing to the secondary body m2
            y-axis is inside the orbital plane of m1 and m2
            z-axis is perpendicular to the orbital plane of m1 and m2
            The x-axis will rotate with the secondary body, making this a non-inertial reference frame.
            This frame is valid only if the two bodies are in circular orbit (will implement for elliptical orbits later)
            This frame is useful for plotting missions from one body to another, as well as to the Lagrange points.
        """
        fig = plt.figure(figsize=(16, 9))
        if plot3d:
            ax = fig.add_subplot(projection='3d')
        else:
            ax = fig.add_subplot()
        
        def draw_orbit(t_clock=0):
            scale = self.rp  # Scale for auxiliary vectors
            if frame == "perifocal":
                h_vec = self.Q_bc_peri(self.t1_clock) @ self.h_vec_bc
                e_vec = self.Q_bc_peri(self.t1_clock) @ self.e_vec_bc
                n_vec = self.Q_bc_peri(self.t1_clock) @ self.n_vec_bc
            elif frame == "perifocal_t0":
                h_vec = self.Q_bc_perit0 @ self.h_vec_bc
                e_vec = self.Q_bc_perit0 @ self.e_vec_bc
                n_vec = self.Q_bc_perit0 @ self.n_vec_bc
            
            ### Scaling vectors for plotting ###
            if self.h == 0:
                h_vec = [0, 0, 0]
            else:
                h_vec = self.h_vec_bc/self.h * scale
            if self.n == 0:
                n_vec = [0, 0, 0]
            else:
                n_vec = self.n_vec_bc/self.n * scale
            if self.e == 0:
                e_vec = [0, 0, 0]
            else:
                e_vec = self.e_vec_bc/self.e * scale

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
                ax.quiver(0, 0, 0, h_vec[0], h_vec[1], h_vec[2], color='g', label='h - Angular Momentum')
                ax.quiver(0, 0, 0, e_vec[0], e_vec[1], e_vec[2], color='r', label='e - Eccentricity')
                ax.quiver(0, 0, 0, n_vec[0], n_vec[1], n_vec[2], color='y', label='n - Ascending Node')
            else:
                ax.plot(r_vecs[:,0], r_vecs[:,1], 'b-', label='Orbit', linewidth=2, zorder=1)
                ax.quiver(0, 0, h_vec[0], h_vec[1], color='g', label='h', scale_units='xy', scale=1, width=0.003)
                ax.quiver(0, 0, e_vec[0], e_vec[1], color='r', label='e', scale_units='xy', scale=1, width=0.003)
                ax.quiver(0, 0, n_vec[0], n_vec[1], color='y', label='n', scale_units='xy', scale=1, width=0.003)

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

            if radius is not None:
                if index == 1:
                    #equator line
                    phi = np.linspace(0, 2 * np.pi, 36)
                    #equator vector in bodycentric frame
                    equator_vec_bc = np.array([x_center + radius * np.cos(phi),
                                            y_center + radius * np.sin(phi),
                                            z_center * np.ones(np.size(phi))])
                    if frame == "perifocal":
                        equator_vec = self.Q_bc_peri(self.t1_clock) @ equator_vec_bc
                    elif frame == "perifocal_t0":
                        equator_vec = self.Q_bc_perit0 @ equator_vec_bc
                    else:
                        equator_vec = equator_vec_bc
                if plot3d:
                    u = np.linspace(0, 2 * np.pi, 36)
                    v = np.linspace(0, np.pi, 36)
                    x = x_center + radius * np.outer(np.cos(u), np.sin(v))
                    y = y_center + radius * np.outer(np.sin(u), np.sin(v))
                    z = z_center + radius * np.outer(np.ones(np.size(u)), np.cos(v))
                    ax.plot_surface(x, y, z, color=color, alpha=0.4, label=label)
                    if index == 1:
                        ax.plot(equator_vec[0], equator_vec[1], equator_vec[2], color='gray', label='Equator')
                else:
                    circle = plt.Circle((x_center, y_center), radius, color=color, alpha=0.4, label=label, zorder=2)
                    plt.gca().add_patch(circle)
                    if index == 1:
                        ax.plot(equator_vec[0], equator_vec[1], color='gray', label='Equator')
            else:
                if plot3d:
                    ax.plot(x_center, y_center, z_center, 'o', color=color, label=label)
                else:
                    plt.plot(x_center, y_center, 'o', color=color, label=label, zorder=2)
            
        def draw_points():
            if self.points_in_orbital_plane:
                for i, point in enumerate(self.points_in_orbital_plane):
                    label = point['name']
                    # List of colors
                    colors = ['red', 'purple', 'brown', 'deeppink', 'olive', 'deepskyblue', 'green']
                    color = colors[i % len(colors)]
                    r_vec = self.r_vec_at_theta(point['theta'], frame=frame)
                    if plot3d:
                        ax.plot(r_vec[0], r_vec[1], r_vec[2], 'o', color=color, label=label)
                    else:
                        ax.plot(r_vec[0], r_vec[1], 'o', color=color, label=label, zorder=2)
        
        def draw_positions():
            if self.positions:
                for i, position in enumerate(self.positions):
                    label = position.get('name', f'Position {i+1}')
                    r_vec, v_vec = self.state_vectors_at_theta(position['theta'], frame=frame)
                    if groundtrack:
                        r_vec = r_vec / np.linalg.norm(r_vec) * self.body1radius
                    if position['t_clock'] is not None:
                        label += f' (t={position['t_clock']:.2f}s)'
                    t_norm = (position['t_clock'] - self.t0_clock) / (self.t1_clock - self.t0_clock)
                    color = self.color_gradient(t_norm)
                    # Plot the secondary body or a marker if no secondary body is provided
                    draw_body(r_vec[0], r_vec[1], r_vec[2], index=2, radius=self.body2radius, color=color, label=label)
                    
                    # Plot velocity vector
                    if velocities:
                        v_vec = v_vec * self.p/20  # Velocity vector scale
                        if groundtrack:
                            v_r = np.dot(v_vec, r_vec/np.linalg.norm(r_vec))
                            v_vec = v_vec - v_r * r_vec/np.linalg.norm(r_vec)
                        if plot3d:
                            ax.quiver(r_vec[0], r_vec[1], r_vec[2], v_vec[0], v_vec[1], v_vec[2], color=color)
                        else:
                            ax.quiver(r_vec[0], r_vec[1], v_vec[0], v_vec[1], color=color, scale_units='xy', scale=1, width=0.003, zorder=5)

        def draw_trajectory():
            if self.trajectory_points:
                annotations = []  # List to store all annotations
                # Trajectory points for annotations
                for point in self.trajectory_points:
                    r_vec = self.r_vec_at_theta(point['theta'], frame=frame)
                    if point['t_clock'] is not None:
                        t_clock = point['t_clock']
                        t_norm = (t_clock - self.t0_clock) / (self.t1_clock - self.t0_clock)
                        color = self.color_gradient(t_norm)
                        if groundtrack:
                            r_vec = r_vec / np.linalg.norm(r_vec) * self.body1radius
                        if plot3d:
                            point_plot = ax.plot(r_vec[0], r_vec[1], r_vec[2], 'o', color=color,
                                                picker=True, pickradius=5)[0]
                        else:
                            point_plot = plt.plot(r_vec[0], r_vec[1], 'o', color=color,
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
                
                # Generate high resolution points for the continuous line
                times = np.linspace(self.t0_clock, self.t1_clock, int((self.t1_clock - self.t0_clock)/60))
                r_vecs = []
                colors = []

                for t in times:
                    r_vec = self.state_vectors_at_t_clock(t, frame=frame)[0]
                    if groundtrack:
                        r_vec = r_vec / np.linalg.norm(r_vec) * self.body1radius
                    t_norm = (t - self.t0_clock) / (self.t1_clock - self.t0_clock)
                    r_vecs.append(r_vec)
                    colors.append(self.color_gradient(t_norm))

                # Convert to numpy array for better performance
                r_vecs = np.array(r_vecs)
                for i in range(len(r_vecs)-1):
                    if plot3d:
                        ax.plot(r_vecs[i:i+2,0], r_vecs[i:i+2,1], r_vecs[i:i+2,2], 
                                '-', color=colors[i], alpha=0.6, linewidth=2)
                    else:
                        plt.plot(r_vecs[i:i+2,0], r_vecs[i:i+2,1], 
                                 '-', color=colors[i], alpha=0.6, linewidth=2, zorder=3)

        if (frame == "perifocal" or
            frame == "perifocal_t0" or
            frame == "bodycentric" or
            frame == "rotatingBodycentric"):
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
        else:
            raise ValueError(f"Invalid frame: {frame}")

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

    def plot_groundtrack(self, frame="rotatingBodycentric"):
        """
        Plots the groundtrack of the orbit.
        """
        fig = plt.figure(figsize=(16, 9))
        ax = fig.add_subplot()

        ########### Plot trajectory points ###########
        if not self.trajectory_points:
            raise ValueError("No trajectory points to plot.")
        
        annotations = []  # List to store all annotations        
        # Trajectory points for annotations
        for point in self.trajectory_points:
            r_vec = self.r_vec_at_theta(point['theta'], frame=frame)
            ra, dec = self.convert_cartesian_to_ra_dec(r_vec)
            if point['t_clock'] is not None:
                t_clock = point['t_clock']
                t_norm = (t_clock - self.t0_clock) / (self.t1_clock - self.t0_clock)
                color = self.color_gradient(t_norm)
                
                point_plot = plt.plot(ra, dec, 'o', color=color,
                                        picker=True, pickradius=5,
                                        markersize=5, zorder=4)[0]
                
                # Create annotation
                annotation = plt.annotate(f't = {t_clock:.2f}s',
                    xy=(ra, dec), xytext=(10, 10),
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
        
        ##### Continuous line #####
        times = np.linspace(self.t0_clock, self.t1_clock, int((self.t1_clock - self.t0_clock)/20))
        for t in times:
            r_vec = self.state_vectors_at_t_clock(t, frame=frame)[0]
            ra, dec = self.convert_cartesian_to_ra_dec(r_vec)
            t_norm = (t - self.t0_clock) / (self.t1_clock - self.t0_clock)
            color = self.color_gradient(t_norm)
            plt.plot(ra, dec, '.', color=color,
                     alpha=0.6, markersize=3, zorder=3)

        
        ########### Plot positions ###########
        if self.positions:
            for i, position in enumerate(self.positions):
                label = position.get('name', f'Position {i+1}')           
                r_vec = self.r_vec_at_theta(position['theta'], frame=frame)
                if position['t_clock'] is not None:
                    label += f' (t={position['t_clock']:.2f}s)'
                t_norm = (position['t_clock'] - self.t0_clock) / (self.t1_clock - self.t0_clock)
                color = self.color_gradient(t_norm)
                ra, dec = self.convert_cartesian_to_ra_dec(r_vec)
                plt.plot(ra, dec, 'o', color=color, label=label,
                         markersize=8, zorder=5)

        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        ax.set_xlabel('Right Ascension (°)')
        ax.set_ylabel('Declination (°)')
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.grid(True)
        ax.axhline(y=0, color='k', linestyle='-', linewidth=1)
        ax.axvline(x=0, color='k', linestyle='-', linewidth=1)

        # Add text for coordinates
        coord_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, 
                            bbox=dict(facecolor='white', alpha=0.7),
                            verticalalignment='top')

        # Show coordinates on hover
        def update_coords(event):
            if event.inaxes == ax:
                ra = event.xdata
                dec = event.ydata
                coord_text.set_text(f'RA: {ra:.2f}°\nDec: {dec:.2f}°')
                plt.draw()

        plt.gcf().canvas.mpl_connect('motion_notify_event', update_coords)

        plt.title(f'Groundtrack')
        plt.tight_layout()
        plt.get_current_fig_manager().window.showMaximized()
        plt.show()

    ####### String representation #######
    def __str__(self):
        """
        Returns a string representation of the orbital parameters.
        """
        # Get all attributes of the class
        attributes = vars(self)
        # Initial string
        result = "Orbital Parameters:\n"
        # Add each attribute to the string
        for name, value in attributes.items():
            # Skip attributes that are None or that should not be shown
            if value is None or name in ['points_in_orbital_plane', 'positions', 'trajectory_points']:
                continue
                
            # Format numbers as float with 3 decimal places
            if isinstance(value, (int, float)):
                result += f"{name}: {value:.3f}\n"
            else:
                result += f"{name}: {value}\n"

        return result




class Orbit(OrbitalBase):
    def __init__(self, mu=None, m1=None, m2=0, a=None, e=None, rp=None, ra=None, h=None, Omega=None, i=None, omega=None, theta0=0, t0_clock=0, body1radius=None, body2radius=None):
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
        """
        super().__init__()
        self.body1radius = body1radius
        self.body2radius = body2radius

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

        self.Omega = Omega
        self.i = i
        self.omega = omega
        self.theta0 = theta0
        self.t0_clock = t0_clock
        self.t1_clock = None
        self.Omega_dot = None
        self.omega_dot = None
        self.omega_body = None
        self.add_zero_state_from_theta(theta0, t0_clock, add_position=True)
        if Omega is not None and i is not None and omega is not None and theta0 is not None:
            self.Q_bc_perit0 = self.Q_from_Euler_angles(Omega, i, omega, pattern="classical")
            self.r_vec_bc, self.v_vec_bc = self.state_vectors_at_theta(theta0, frame="bodycentric")
            self.h_vec_bc = self.h_vec_from_state_vectors(self.r_vec_bc, self.v_vec_bc)
            self.e_vec_bc = self.e_vec_from_mu_state_vectors(self.mu, self.r_vec_bc, self.v_vec_bc)
            self.n_vec_bc = self.n_vec_from_h_vec(self.h_vec_bc)
            self.r = np.linalg.norm(self.r_vec_bc)
            self.v = np.linalg.norm(self.v_vec_bc)
            self.n = np.linalg.norm(self.n_vec_bc)

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
    def init_from_state_vectors(cls, r_vec_bc, v_vec_bc, mu=None, m1=None, m2=0, t0_clock=0, body1radius=None, body2radius=None):
        """
        Given a position and velocity vectors represented in some inertial frame of reference
        centered in the primary body, creates a new instance of Orbit.

        :param r_vec_bc: Position vector in the bodycentric frame of reference (km)
        :param v_vec_bc: Velocity vector in the bodycentric frame of reference (km/s)
        :param mu: Standard gravitational parameter (km³/s²)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param t0_clock: Time of the starting point (s)
        :param body1radius: Radius of the central body (km) (optional)
        :param body2radius: Radius of the orbiting body (km) (optional)
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Provide mu or m1 (and optionally m2) to calculate mu.")
        
        # Calculate all orbital parameters
        orbital_parameters = cls.orbital_parameters_from_state_vectors(r_vec_bc, v_vec_bc, mu)
        h = orbital_parameters['h']
        e = orbital_parameters['e']
        Omega = orbital_parameters['Omega']
        i = orbital_parameters['i']
        omega = orbital_parameters['omega']
        theta0 = orbital_parameters['theta0']
        
        return cls(mu=mu, m1=m1, m2=m2, e=e, h=h, Omega=Omega, i=i, omega=omega, theta0=theta0, t0_clock=t0_clock, body1radius=body1radius, body2radius=body2radius)

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
    