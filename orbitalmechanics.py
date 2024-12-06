import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar, fsolve

class OrbitalBase:
    """
    Base class for all orbital mechanics classes.
    """
    G = 6.67430e-20  # Constante gravitacional em km^3/kg/s^2

    #####   mu   #####
    def _calc_mu(self, m1, m2=0):
        """
        Calculates the gravitational parameter mu based on the masses of two bodies.
        
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body, default is 0 (kg)
        """
        return self.G * (m1 + m2)

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
        self.body1radius = body1radius
        self.body2radius = body2radius

        if mu is not None:
            self.mu = mu  # Gravitational parameter
        elif m1 is not None:
            self.mu = self._calc_mu(m1, m2)
        else:
            raise ValueError("Provide either mu or m1 (and optionally m2) to calculate mu.")
        
        if rp is not None and e is not None:
            self._rp = rp
            self._e = e
            self._h = self._calc_h_from_mu_rp_e(self.mu, self._rp, self._e)
            self._ra = self._calc_ra_from_mu_h_e(self.mu, self._h, self._e)
            self._a = self._calc_a_from_rp_ra(self._rp, self._ra)
        elif rp is not None and ra is not None:
            self._rp = rp
            self._ra = ra
            self._a = self._calc_a_from_rp_ra(self._rp, self._ra)
            self._e = self._calc_e_from_rp_ra(self._rp, self._ra)
            self._h = self._calc_h_from_mu_rp_e(self.mu, self._rp, self._e)
        elif rp is not None and h is not None:
            self._rp = rp
            self._h = h
            self._e = self._calc_e_from_mu_h_rp(self.mu, self._h, self._rp)
            self._ra = self._calc_ra_from_mu_h_e(self.mu, self._h, self._e)
            self._a = self._calc_a_from_rp_ra(self._rp, self._ra)
        elif e is not None and h is not None:
            self._e = e
            self._h = h
            self._rp = self._calc_rp_from_mu_h_e(self.mu, self._h, self._e)
            self._ra = self._calc_ra_from_mu_h_e(self.mu, self._h, self._e)
            self._a = self._calc_a_from_rp_ra(self._rp, self._ra)
        elif a is not None and e is not None:
            if e == 1:
                raise ValueError("Can't find rp from 'a' and 'e' for parabolic (e=1) orbits.")
            self._a = a
            self._e = e
            self._rp = self._calc_rp_from_a_e(self._a, self._e)
            self._ra = self._calc_ra_from_a_e(self._a, self._e)
            self._h = self._calc_h_from_mu_rp_e(self.mu, self._rp, self._e)
        else:
            raise ValueError("Provide sufficient parameters to define the orbit (rp and e, rp and ra, rp and h, e and h or a and e).")
        
        self.b = self._calc_b_from_a_e(self._a, self._e)
        self.epsilon = self._calc_epsilon_from_mu_h_e(self.mu, self._h, self._e)
        self.T = self._calc_T_from_mu_a(self.mu, self._a)
        self.p = self._calc_p_from_mu_h(self.mu, self._h)
    
    @classmethod
    def init_from_2points(cls, mu=None, m1=None, m2=0, r1=None, theta1=None, r2=None, theta2=None, body1radius=None, body2radius=None):
        """
        Creates a new instance of Orbit from two points (r,θ).
        
        :param mu: Standard gravitational parameter (km³/s²)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param r1: Radius of the first point (km)
        :param theta1: True anomaly of the first point (degrees)
        :param r2: Radius of the second point (km)
        :param theta2: True anomaly of the second point (degrees)
        :return: New instance of Orbit
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Forneça mu ou m1 (e opcionalmente m2) para calcular mu.")
        
        # Converter ângulos para radianos
        theta1 = np.radians(theta1)
        theta2 = np.radians(theta2)
        
        # Calcular e
        e = (r2/r1 - 1) / (np.cos(theta1) - (r2/r1)*np.cos(theta2))
        
        # Calcular h
        h = np.sqrt(mu * r1 * (1 + e * np.cos(theta1)))
        
        # Simplificado para sempre usar mu
        return cls(mu=mu, e=e, h=h, body1radius=body1radius, body2radius=body2radius)

    @classmethod
    def init_from_r_v_gamma(cls, mu=None, m1=None, m2=0, r=None, v=None, gamma=None, body1radius=None, body2radius=None):
        """
        Creates a new instance of Orbit from distance, velocity and launch angle.

        :param mu: Standard gravitational parameter (km³/s²)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param r: Distance from the primary body center to the point (km)
        :param v: Velocity at the point (km/s)
        :param gamma: Launch angle (degrees)
        :return: New instance of Orbit
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Forneça mu ou m1 (e opcionalmente m2) para calcular mu.")
        
        gamma = np.radians(gamma)
        vr = v * np.sin(gamma)
        vt = v * np.cos(gamma)
        h = r * vt
        e = np.sqrt(vr**2*h**2/mu**2 + (h**2/(mu*r) - 1)**2)
        return cls(mu=mu, e=e, h=h, body1radius=body1radius, body2radius=body2radius)

    @classmethod
    def init_from_perifocal_vectors(cls, mu=None, m1=None, m2=0, r_vec=None, v_vec=None, body1radius=None, body2radius=None):
        """
        Creates a new instance of Orbit from position and velocity vectors in the perifocal frame of reference.

        :param mu: Standard gravitational parameter (km³/s²)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param r_vec: Position vector in the perifocal frame of reference (km)
        :param v_vec: Velocity vector in the perifocal frame of reference (km/s)
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Forneça mu ou m1 (e opcionalmente m2) para calcular mu.")
        
        h_vec = np.cross(r_vec, v_vec)
        h = np.linalg.norm(h_vec)
        theta = np.arctan2(r_vec[1], r_vec[0])
        e = (h**2 / (mu * np.linalg.norm(r_vec)) - 1)/np.cos(theta)
        return cls(mu=mu, e=e, h=h, body1radius=body1radius, body2radius=body2radius)
    
    @classmethod
    def init_from_r_vec_v_vec(cls, mu=None, m1=None, m2=0, r_vec=None, v_vec=None, name=None, body1radius=None, body2radius=None):
        """
        Given a position and velocity vectors represented in some inertial frame of reference
        centered in the primary body, creates a new instance of Orbit.

        :param mu: Standard gravitational parameter (km³/s²)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param r_vec: Position vector in the inertial frame of reference (km)
        :param v_vec: Velocity vector in the inertial frame of reference (km/s)
        :param name: Name of the body
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Forneça mu ou m1 (e opcionalmente m2) para calcular mu.")
        
        r_vec_0 = r_vec
        v_vec_0 = v_vec
        r0 = np.linalg.norm(r_vec_0)
        vr0 = np.dot(r_vec_0, v_vec_0) / r0
        h = np.linalg.norm(np.cross(r_vec_0, v_vec_0))

        #f, g, f_dot, g_dot = cls.lagrange_coefficients(r_vec_0, v_vec_0)
        #delta_theta = 120
        #r = f(delta_theta)*r_vec_0 + g(delta_theta)*v_vec_0
        #v = f_dot(delta_theta)*r_vec_0 + g_dot(delta_theta)*v_vec_0

        e = np.sqrt(vr0**2 * h**2 / mu**2 + (h**2 / (mu * r0) - 1)**2)
        rp = h**2 / (mu * (1 + e))
        
        # Criar a instância da órbita
        orbit_instance = cls(mu=mu, e=e, h=h, rp=rp, body1radius=body1radius, body2radius=body2radius)
        
        # Adicionar o corpo usando o método add_body
        theta = np.degrees(np.arctan2(r_vec[1], r_vec[0]))  # Calcula a anomalia verdadeira
        orbit_instance.add_body(name=name, theta=theta)
        
        return orbit_instance

    @classmethod
    def init_from_r_v_theta(cls, mu=None, m1=None, m2=0, r=None, v=None, theta=None, body1radius=None, body2radius=None):
        """
        Creates a new instance of Orbit from distance, velocity and true anomaly.

        :param mu: Standard gravitational parameter (km³/s²)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param r: Distance from the primary body center to the point (km)
        :param v: Velocity at the point (km/s)
        :param theta: True anomaly (degrees)
        """
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
        if mu is None:
            raise ValueError("Forneça mu ou m1 (e opcionalmente m2) para calcular mu.")
        
        theta = np.radians(theta)

        # solve second degree equation for eccentricity
        u = r*v**2/mu
        c1 = 1
        c2 = np.cos(theta)*(2 - u)
        c3 = 1 - u
        delta = c2**2 - 4*c1*c3
        e = (-c2 + np.sqrt(delta))/(2*c1)
                
        # calculate specific angular momentum
        h = np.sqrt(mu * r * (1 + e*np.cos(theta)))

        return cls(mu=mu, e=e, h=h, body1radius=body1radius, body2radius=body2radius)

    def add_position(self, theta, name=None):
        """
        Adds an orbital position to the orbit.
        """
        r = self.r_at_theta(theta)
        v = self.v_at_r(r)
        self.positions = getattr(self, 'positions', [])
        position = {
            'nome': name,
            'r': r,
            'v': v, 
            'theta': theta
        }
        self.positions.append(position)
        
    ############# Calculations #############
    #####   rp   #####
    def _calc_rp_from_a_e(self, a, e):
        """
        Calculates the periapsis distance using semi-major axis and eccentricity.
        Note: For hyperbolic orbits (e > 1), 'a' must be a negative value.
        """
        return a * (1 - e)

    def _calc_rp_from_mu_h_e(self, mu, h, e):
        """
        Calculates the periapsis distance using gravitational parameter, specific angular momentum and eccentricity.
        """
        return h**2 / (mu * (1 + e))

    #####   ra   #####
    def _calc_ra_from_a_e(self, a, e):
        """
        Calculates the apoapsis distance using semi-major axis and eccentricity.
        Note: For hyperbolic orbits (e > 1), 'a' must be a negative value, 
        which will result in a negative distance.
        """
        if e == 1:
            return float('inf')
        return a * (1 + e)
    
    def _calc_ra_from_mu_h_e(self, mu, h, e):
        """
        Calculates the apoapsis distance using gravitational parameter, specific angular momentum and eccentricity.
        """
        if e == 1:
            return float('inf')
        else:
            return h**2 / (mu * (1 - e))

    #####   a   #####
    def _calc_a_from_rp_ra(self, rp, ra):
        """
        Calculates the semi-major axis using periapsis and apoapsis distances.
        Note: For hyperbolic orbits (e > 1), 'ra' must be a negative value, 
        which will result in a negative distance.
        """
        if ra == float('inf'):
            return float('inf')
        return (rp + ra) / 2

    #####   e   #####
    def _calc_e_from_rp_ra(self, rp, ra):
        """
        Calculates the eccentricity using periapsis and apoapsis distances.
        Note: For hyperbolic orbits (e > 1), 'ra' must be a negative value. 
        """
        if ra == float('inf'):
            return 1
        return (ra - rp) / (ra + rp)
    
    def _calc_e_from_mu_h_rp(self, mu, h, rp):
        """
        Calculates the eccentricity using gravitational parameter, specific angular momentum and periapsis distance.
        """
        return h**2 / (mu * rp) - 1
    
    #####   h   #####
    def _calc_h_from_mu_a_e(self, mu, a, e):
        """
        Calculates the specific angular momentum using gravitational parameter, semi-major axis and eccentricity.
        For hyperbolic orbits (e > 1), a should be negative.
        """
        if e == 1:
            raise ValueError("Can't calculate h from 'a' and 'e' for parabolic (e=1) orbits.")
        return np.sqrt(mu * a * (1 - e**2))
    
    def _calc_h_from_mu_rp_e(self, mu, rp, e):
        """
        Calculates the specific angular momentum using gravitational parameter, periapsis distance and eccentricity.
        """
        return np.sqrt(mu * rp * (1 + e))

    #####   b   #####
    def _calc_b_from_a_e(self, a, e):
        """
        Calculates the semi-minor axis using semi-major axis and eccentricity.
        Note: For hyperbolic orbits (e > 1), 'a' must be a negative value.
        """
        if e == 1:
            return float('inf')
        elif e > 1:
            return -a * np.sqrt(e**2 - 1)
        else:
            return a * np.sqrt(1 - e**2)
    
    #####   p   #####
    def _calc_p_from_mu_h(self, mu, h):
        """
        Calculates the semi-latus rectum using gravitational parameter and specific angular momentum.
        """
        return h**2 / mu

    #####   epsilon   #####
    def _calc_epsilon_from_mu_h_e(self, mu, h, e):
        """
        Calculates the specific orbital energy using gravitational parameter,
        specific angular momentum and eccentricity.
        """
        return -mu**2 / (2 * h**2) * (1 - e**2)
    
    def _calc_epsilon_from_mu_a(self, mu, a):
        """
        Calculates the specific orbital energy using gravitational parameter and semi-major axis.
        Note: For hyperbolic orbits (e > 1), 'a' must be a negative value.
        """
        return -mu / (2 * a)

    #####   T   #####
    def _calc_T_from_mu_a(self, mu, a):
        """
        Calculates the orbital period using gravitational parameter and semi-major axis.
        Note: For hyperbolic orbits (e > 1), 'a' must be a negative value.
        """
        if a == float('inf') or a < 0:
            return float('inf')
        return 2 * np.pi * np.sqrt(a**3 / mu)

    def _update(self, param):
        """
        Recalculates the orbital parameters based on the changed parameter.
        TODO: test for each parameter
        """
        if param in ['rp', 'ra']:
            self._a = self._calc_a_from_rp_ra(self._rp, self._ra)
            self._e = self._calc_e_from_rp_ra(self._rp, self._ra)
            
        if param in ['a', 'e']:
            self._rp = self._calc_rp_from_a_e(self._a, self._e)
            self._ra = self._calc_ra_from_a_e(self._a, self._e)
            
        self._b = self._calc_b_from_a_e(self._a, self._e)
        self.epsilon = self._calc_epsilon_from_a(self._a)
        self._h = self._calc_h_from_a_e(self._a, self._e)
        self.T = self._calc_T_from_mu_a(self.mu, self._a)
        self.r_mean = self._calc_r_mean_from_rp_ra(self._rp, self._ra)
    
    def __str__(self):
        """
        Returns a string representation of the orbital parameters.
        """
        return (f"Orbital Parameters:\n"
                f"Semi-major axis: {self._a:.3f} km\n"
                f"Semi-minor axis: {self.b:.3f} km\n"
                f"Eccentricity: {self._e:.4f} ({self.type()})\n"
                f"Periapsis: {self._rp:.3f} km\n"
                f"Apoapsis: {self._ra:.3f} km\n"
                f"Orbital Period: {self.T:.3f} s\n"
                f"Specific Orbital Energy: {self.epsilon:.3f} km^2/s^2\n"
                f"Specific Angular Momentum: {self._h:.3f} km^2/s\n"
                f"Semi-latus Rectum: {self.p:.3f} km\n")
    
    def type(self):
        """
        Classifies the type of orbit based on eccentricity.
        :return: Type of orbit as a string ("Circular", "Elliptical", "Parabolic", "Hyperbolic")
        """
        if abs(self._e) < 1e-10:  # very close to zero
            return "Circular"
        elif 0 < self._e < 1:
            return "Elliptical"
        elif self._e == 1:
            return "Parabolic"
        elif self._e > 1:
            return "Hyperbolic"
        else:
            return "Invalid"

    def theta_averaged_r(self):
        """
        Calculates the theta-averaged radius using periapsis and apoapsis distances.
        Note: Can't be used for parabolic or hyperbolic orbits.
        """
        if self._e >= 1:
            raise ValueError("Can't calculate theta-averaged radius for parabolic or hyperbolic orbits.")
        return (self._rp * self._ra)**(1/2)
    
    def time_averaged_r(self):
        """
        Calculates the time-averaged radius using semi-major axis and eccentricity.
        Note: Can't be used for parabolic or hyperbolic orbits.
        """
        if self._e >= 1:
            raise ValueError("Can't calculate time-averaged radius for parabolic or hyperbolic orbits.")
        return self._a * (1 + self._e**2/2)

    def v_at_r(self, r):
        """
        Calculates the orbital velocity at a given distance using the vis-viva equation.
        
        :param r: Distance from the primary body center to the point (km)
        :return: Velocity at the point (km/s)
        """
        if self._e == 1:
            return np.sqrt(self.mu * 2 / r)
        else:
            return np.sqrt(self.mu * (2 / r - 1 / self._a))
    
    def r_at_theta(self, theta):
        """
        Calculates the distance from the primary body center to a point on the orbit at a given true anomaly.
        """
        theta = np.radians(theta)
        return self._h**2 / (self.mu) * 1/(1 + self._e * np.cos(theta))
    
    def theta_at_r(self, r):
        """
        Calculates the true anomaly at a given distance from the primary body center.
        Returns both possible values as a list, with angles between 0 and 360 degrees.
        """
        theta1 = np.arccos((self._h**2 / (self.mu * r) - 1) / self._e)
        theta2 = 2*np.pi - theta1
        return [np.degrees(theta1), np.degrees(theta2)]
    
    def gamma_at_theta(self, theta):
        """
        Calculates the flight path angle at a given true anomaly.
        """
        theta = np.radians(theta)
        return np.degrees(np.arctan(self._e * np.sin(theta) / (1 + self._e * np.cos(theta))))
    
    def v_esc(self, r):
        """
        Calculates the escape velocity at a given distance.
        """
        return np.sqrt(2 * self.mu / r)
    
    def v_inf(self):
        """
        Calculates the hyperbolic excess velocity.
        """
        return np.sqrt(self.mu/self._a)
    
    def C3(self):
        """
        Calculates the hyperbolic energy.
        """
        return self.epsilon * 2
    
    def turn_angle(self):
        """
        Calculates the turn angle.
        """
        return np.degrees(2 * np.arcsin(1/self._e))
    
    def theta_inf(self):
        """
        Calculates the true anomaly at infinity.
        """
        return np.degrees(np.arccos(-1/self._e))
    
    def aiming_radius(self):
        """
        Calculates the aiming radius.
        """
        return self.b
    
    def r_vec_perifocal(self, theta):
        """
        Calculates the position vector in the perifocal frame of reference.
        """
        theta = np.radians(theta)
        return self._h**2 / (self.mu * (1 + self._e*np.cos(theta))) * np.array([np.cos(theta), np.sin(theta), 0])
    
    def v_vec_perifocal(self, theta):
        """
        Calculates the velocity vector in the perifocal frame of reference.
        """
        theta = np.radians(theta)
        return self.mu / self._h * np.array([-np.sin(theta), self._e + np.cos(theta), 0])
    
    def state_perifocal(self, theta):
        """
        Calculates the state vector in the perifocal frame of reference.
        """
        r = self.r_vec_perifocal(theta)
        v = self.v_vec_perifocal(theta)
        return r, v

    def E_from_theta(self, theta):
        """
        Calculates the eccentric anomaly using the true anomaly.
        """
        theta = np.radians(theta)
        if self._e < 1:
            return np.degrees(2 * np.arctan(np.sqrt((1 - self._e) / (1 + self._e)) * np.tan(theta/2)))
        elif self._e > 1:
            return np.degrees(2 * np.arctanh(np.sqrt((self._e - 1) / (self._e + 1)) * np.tan(theta/2)))
        else:
            raise ValueError("Can't calculate eccentric anomaly from true anomaly for parabolic orbits.")

    def M_from_E(self, E):
        """
        Calculates the mean anomaly using the eccentric anomaly.
        """
        E = np.radians(E)
        if self._e < 1:
            return np.degrees(E - self._e * np.sin(E))
        elif self._e > 1:
            return np.degrees(self._e * np.sinh(E) - E)
        else:
            raise ValueError("Can't calculate mean anomaly from eccentric anomaly for parabolic orbits.")
    
    def M_from_theta(self, theta):
        """
        Calculates the mean anomaly using the true anomaly.
        """
        if self._e > 1 or self._e < 1:
            raise ValueError("Can't calculate mean anomaly from true anomaly for elliptical or hyperbolic orbits.")
        
        theta = np.radians(theta)
        return np.degrees(1/6 * np.tan(theta/2)**3 + 1/2 * np.tan(theta/2))
    
    def t_from_M(self, M):
        """
        Calculates the time using the mean anomaly.
        """
        M = np.radians(M)
        if self._e < 1:
            return M * self.T / (2 * np.pi)
        elif self._e == 1:
            return M * self._h**3 / self.mu**2
        else:
            return M * self._h**3 / (self.mu**2 * (self._e**2 - 1)**(3/2))
    
    def t_from_theta(self, theta):
        """
        Calculates the time using the true anomaly.
        """
        if self._e < 1 or self._e > 1: # elliptical or hyperbolic
            E = self.E_from_theta(theta)
            M = self.M_from_E(E)
            return self.t_from_M(M)
        else: # parabolic
            M = self.M_from_theta(theta)
            return self.t_from_M(M)

    def M_from_t(self, t):
        """
        Calculates the mean anomaly using the time.
        """
        if self._e < 1:
            return np.degrees(2 * np.pi * t / self.T)
        elif self._e == 1:
            return np.degrees(self.mu**2 * t/ self._h**3)
        else:
            return np.degrees(self.mu**2 / self._h**3 * (self._e**2 - 1)**(3/2) * t)

    def E_from_M(self, M):
        """
        Calculates the eccentric anomaly using the mean anomaly.
        """
        if self._e == 1:
            raise ValueError("Can't calculate eccentric anomaly for parabolic orbits.")
        
        M = np.radians(M)
        if self._e < 1:
            if M < np.pi:
                E0 = M + self._e/2
            else:
                E0 = M - self._e/2
        else:
            E0 = M

        def f(E, M):
            if self._e < 1:
                return E - self._e * np.sin(E) - M
            else:
                return self._e * np.sinh(E) - E - M

        def f_prime(E, M):
            if self._e < 1:
                return 1 - self._e * np.cos(E)
            else:
                return self._e * np.cosh(E) - 1
        
        E = root_scalar(f, fprime=f_prime, x0=E0, args=(M), method='newton').root

        return np.degrees(E)

    def theta_from_E(self, E):
        """
        Calculates the true anomaly using the eccentric anomaly.
        """
        E = np.radians(E)
        if self._e < 1:
            return np.degrees(2 * np.arctan(np.sqrt((1 + self._e) / (1 - self._e)) * np.tan(E/2)))
        elif self._e > 1:
            return np.degrees(2 * np.arctan(np.sqrt((self._e + 1) / (self._e - 1)) * np.tanh(E/2)))
        else:
            raise ValueError("Can't calculate true anomaly from eccentric anomaly for parabolic orbits.")

    def theta_from_t(self, t):
        """
        Calculates the true anomaly using the time.
        """
        if self._e < 1 or self._e > 1: # elliptical or hyperbolic
            M = self.M_from_t(t)
            E = self.E_from_M(M)
            return np.mod(self.theta_from_E(E), 360)
        else: # parabolic
            M = self.M_from_t(t)
            M = np.radians(M)
            z = (3*M + np.sqrt(1 + 9*M**2))**(1/3)
            theta = np.degrees(2 * np.arctan(z - 1/z))
            return np.mod(theta, 360)

    def lagrange_coefficients(self, r_vec_0, v_vec_0):
        """
        Calcula e retorna as funções dos coeficientes de Lagrange.
        """
        r0 = np.linalg.norm(r_vec_0)
        v0 = np.linalg.norm(v_vec_0)
        h = np.linalg.norm(np.cross(r_vec_0, v_vec_0))
        vr0 = np.dot(r_vec_0, v_vec_0) / r0

        def r_func(delta_theta):
            delta_theta = np.radians(delta_theta)
            return h**2 / (self.mu * (1 + (h**2/(self.mu*r0) - 1)*np.cos(delta_theta) - h*vr0/self.mu * np.sin(delta_theta)))

        def f(delta_theta):
            delta_theta = np.radians(delta_theta)
            r = r_func(delta_theta)
            return 1 - self.mu*r/h**2 * (1 - np.cos(delta_theta))

        def g(delta_theta):
            delta_theta = np.radians(delta_theta)
            r = r_func(delta_theta)
            return r*r0/h * np.sin(delta_theta)

        def f_dot(delta_theta):
            delta_theta = np.radians(delta_theta)
            r = r_func(delta_theta)
            return self.mu/h * (1 - np.cos(delta_theta))/np.sin(delta_theta) * (self.mu/h**2 * (1 - np.cos(delta_theta)) - 1/r0 - 1/r)

        def g_dot(delta_theta):
            delta_theta = np.radians(delta_theta)
            r = r_func(delta_theta)
            return 1 - self.mu*r0/h**2 * (1 - np.cos(delta_theta))

        return f, g, f_dot, g_dot

    @property
    def a(self):
        return self._a
    
    @a.setter 
    def a(self, value):
        self._a = value
        self._update('a')
        
    @property
    def e(self):
        return self._e
    
    @e.setter
    def e(self, value):
        self._e = value
        self._update('e')
        
    @property
    def rp(self):
        return self._rp
    
    @rp.setter
    def rp(self, value):
        self._rp = value
        self._update('rp')
        
    @property
    def ra(self):
        return self._ra
    
    @ra.setter
    def ra(self, value):
        self._ra = value
        self._update('ra')

    @property
    def h(self):
        return self._h
    
    @h.setter
    def h(self, value):
        self._h = value
        self._update('h')

    def plot(self, points=None, plot_positions=False):
        """
        Plots the orbit with optional points.

        :param points: List of points to plot (r, theta)
        :param plot_positions: Boolean to plot the already added orbital positions
        """
        # Criar array de anomalia verdadeira
        if self._e == 1:
            theta_min = -120
            theta_max = 120
            if plot_positions:
                for pos in self.positions:
                    if pos['theta'] < theta_min:
                        theta_min = pos['theta']
                    if pos['theta'] > theta_max:
                        theta_max = pos['theta']
            theta = np.linspace(theta_min, theta_max, 1000)
        elif self._e > 1:
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
        plt.plot(x, y, 'b-', label='Órbita')
        
        # Plotar o corpo central com o raio correto se fornecido
        if self.body1radius is not None:
            circle = plt.Circle((0, 0), self.body1radius, color='r', alpha=0.3, label='Central Body')
            plt.gca().add_patch(circle)
        else:
            plt.plot(0, 0, 'ro', label='Corpo Central')
            
        # Plotar pontos adicionais se fornecidos
        if points is not None:
            for point in points:
                r, theta = point
                theta = np.radians(theta)
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                plt.plot(x, y, 'go')
        
        if plot_positions:
            # Lista de cores para os pontos
            cores = ['red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'blue']
            for i, position in enumerate(self.positions):
                theta = np.radians(position['theta'])
                x = position['r'] * np.cos(theta)
                y = position['r'] * np.sin(theta)
                # Plotar o corpo 2 com o raio correto se fornecido
                label = position.get('nome', f'Ponto {i+1}')
                cor = cores[i % len(cores)]  # Usa o módulo para reciclar cores se houver mais pontos que cores
                if self.body2radius is not None:
                    circle2 = plt.Circle((x, y), self.body2radius, color=cor, alpha=0.3, label=label)
                    plt.gca().add_patch(circle2)
                else:
                    plt.plot(x, y, 'o', color=cor, label=label)

        plt.grid(True)
        plt.axis('equal')
        plt.xlabel('x (km)')
        plt.ylabel('y (km)')
        plt.title('Órbita')
        plt.legend()
    
    ############# Geocentric Equatorial Frame #############
    @classmethod
    def ra_dec_from_r_vec(cls, r_vec):
        """
        Calculates the right ascension and declination from the position vector.
        """
        r = np.linalg.norm(r_vec)
        ra = np.mod(np.degrees(np.arctan2(r_vec[1], r_vec[0])), 360)
        dec = np.mod(np.degrees(np.arcsin(r_vec[2]/r)), 360)
        return ra, dec

class Lagrange_mechanics(OrbitalBase):
    """
    Class to calculate the position and velocity vectors after a given true anomaly change using Lagrange coefficients.
    """
    def __init__(self, mu=None, m1=None, m2=0, r_vec_0=None, v_vec_0=None, delta_theta=None, body1radius=None, body2radius=None):
        if mu is not None:
            self.mu = mu  # Gravitational parameter
        elif m1 is not None:
            self.mu = self._calc_mu(m1, m2)
        else:
            raise ValueError("Provide either mu or m1 (and optionally m2) to calculate mu.")
        
        self.body1radius = body1radius
        self.body2radius = body2radius

        if r_vec_0 is not None and v_vec_0 is not None and delta_theta is not None:
            self.r_vec_0 = r_vec_0
            self.v_vec_0 = v_vec_0
            self.delta_theta = delta_theta
        else:
            raise ValueError("Provide r_vec_0, v_vec_0 and delta_theta.")
        
        self.f, self.g, self.f_dot, self.g_dot = self.calc_lagrange_coefficients(self.r_vec_0, self.v_vec_0, self.delta_theta)
        self.r_vec, self.v_vec = self.r_vec_v_vec_after_delta_theta(self.delta_theta)
        self._h = np.linalg.norm(np.cross(self.r_vec_0, self.v_vec_0))
        r0 = np.linalg.norm(self.r_vec_0)
        vr0 = np.dot(self.r_vec_0, self.v_vec_0) / r0
        self._e = np.sqrt(vr0**2 * self._h**2 / self.mu**2 + (self._h**2 / (self.mu * r0) - 1)**2)
        self._rp = self._h**2 / (self.mu * (1 + self._e))
        if self._e == 1:
            self._ra = float('inf')
        else:
            self._ra = self._h**2 / (self.mu * (1 - self._e))
        
        if self._ra == float('inf'):
            self._a = float('inf')
        else:
            self._a = (self._rp + self._ra) / 2
        
        if self._e > 1:
            self.theta_inf = np.degrees(np.arccos(-1/self._e))
        
        if vr0 < 0:
            self.theta0 = self.theta_at_r(r0)[1]
        else:
            self.theta0 = self.theta_at_r(r0)[0]
        self.r0 = r0
        r = np.linalg.norm(self.r_vec)
        vr = np.dot(self.r_vec, self.v_vec) / r
        if vr < 0:
            self.theta = self.theta_at_r(r)[1]
        else:
            self.theta = self.theta_at_r(r)[0]
        self.r = r

    ############# Calculations #############
    #####   Lagrange coefficients   #####
    def calc_lagrange_coefficients(self, r_vec_0, v_vec_0, delta_theta):
        """
        Calculates the Lagrange coefficients.
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
    
    def r_vec_v_vec_after_delta_theta(self, delta_theta):
        """
        Calculates the position and velocity vectors after a given true anomaly change.
        """
        r_vec = self.f*self.r_vec_0 + self.g*self.v_vec_0
        v_vec = self.f_dot*self.r_vec_0 + self.g_dot*self.v_vec_0
        return r_vec, v_vec

    def theta_at_r(self, r):
        """
        Calculates the true anomaly at a given distance from the primary body center.
        Returns both possible values as a list, with angles between 0 and 360 degrees.
        """
        theta1 = np.arccos((self._h**2 / (self.mu * r) - 1) / self._e)
        theta2 = 2*np.pi - theta1
        return [np.degrees(theta1), np.degrees(theta2)]

    def r_at_theta(self, theta):
        """
        Calculates the distance from the primary body center to a point on the orbit at a given true anomaly.
        """
        theta = np.radians(theta)
        return self._h**2 / (self.mu) * 1/(1 + self._e * np.cos(theta))

    def plot(self, points=None):
        """
        Plots the orbit with optional points.
        """
        # Criar array de anomalia verdadeira
        if self._e == 1:
            theta = np.linspace(-120, 120, 1000)
        elif self._e > 1:
            #theta = np.linspace(np.radians(-self.theta_inf()), np.radians(self.theta_inf()), 1000)
            epsilon = 15
            theta = np.linspace(-self.theta_inf + epsilon, self.theta_inf - epsilon, 1000)
        else:
            theta = np.linspace(0, 360, 1000)
        
        # Calcular raio para cada ângulo
        r = self.r_at_theta(theta)
        
        # Converter para coordenadas cartesianas
        x = r * np.cos(np.radians(theta))
        y = r * np.sin(np.radians(theta))
        
        # Criar o plot
        plt.figure(figsize=(8, 8))
        plt.plot(x, y, 'b-', label='Órbita')
        
        # Plotar o corpo central com o raio correto se fornecido
        if self.body1radius is not None:
            circle = plt.Circle((0, 0), self.body1radius, color='r', alpha=0.3, label='Central Body')
            plt.gca().add_patch(circle)
        else:
            plt.plot(0, 0, 'ro', label='Corpo Central')
            
        # Plotar o corpo secundário com o raio correto se fornecido
        if self.body2radius is not None:
            circle2 = plt.Circle((self.a, 0), self.body2radius, color='g', alpha=0.3, label='Orbiting Body')
            plt.gca().add_patch(circle2)
        
        # Plotar pontos adicionais se fornecidos
        if points is not None:
            for point in points:
                r, theta = point
                x = r * np.cos(np.radians(theta))
                y = r * np.sin(np.radians(theta))
                plt.plot(x, y, 'go')
        
        plt.grid(True)
        plt.axis('equal')
        plt.xlabel('x (km)')
        plt.ylabel('y (km)')
        plt.title('Órbita')
        plt.legend()

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
        Calculates the velocity vector from some specific Jacobi constant.
        """
        r1 = np.linalg.norm(np.array([r_vec[0]+self.pi2*self.r12, r_vec[1], r_vec[2]]))
        r2 = np.linalg.norm(np.array([r_vec[0]-self.pi1*self.r12, r_vec[1], r_vec[2]]))
        v = np.sqrt(2*C + self.omega**2 * (r_vec[0]**2 + r_vec[1]**2) + 2*self.mu1/r1 + 2*self.mu2/r2)
        return v

    def trajectory(self, r_vec_0, v_vec_0, t_span, t_eval=None, method='RK45', plot=False):
        """
        Simulates the trajectory of the spacecraft.
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

    def plot(self, points):
            """
            Plots the Lagrange points.
            """ 
            plt.figure(figsize=(10, 10))
            
            if self.body1radius is not None:
                circle = plt.Circle((-self.pi2*self.r12, 0), self.body1radius, color='b', alpha=0.3, label='Corpo Primário')
                plt.gca().add_patch(circle)
            else:
                plt.plot(-self.pi2*self.r12, 0, 'bo', label='Corpo Primário')
            if self.body2radius is not None:
                circle = plt.Circle(((1-self.pi2)*self.r12, 0), self.body2radius, color='r', alpha=0.3, label='Corpo Secundário')
                plt.gca().add_patch(circle)
            else:
                plt.plot((1-self.pi2)*self.r12, 0, 'ro', label='Corpo Secundário')
            
            for point in points:
                plt.plot(point[0], point[1], 'k.', markersize=3)
            
            plt.grid(True)
            plt.axis('equal')
            plt.xlabel('x (km)')
            plt.ylabel('y (km)') 
            plt.legend()
            plt.show()
        
class Trajectory(OrbitalBase):
    """
    Class to calculate the trajectory of a spacecraft.
    """
    def __init__(self, r_vec_0, v_vec_0, t0=0, mu=None, m1=None, m2=0, body1radius=None, body2radius=None):
        """
        Initializes the Trajectory class.
        
        In Body-Centric Equatorial Frame.
        :param mu: Gravitational parameter (km^3/s^2)
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body (kg), default is 0
        :param r_vec_0: Initial position vector (km)
        :param v_vec_0: Initial velocity vector (km/s)
        :param t0: Initial time (s), default is 0
        :param body1radius: Radius of the primary body (km), default is None
        :param body2radius: Radius of the secondary body (km), default is None
        """
        self.body1radius = body1radius
        self.body2radius = body2radius

        if mu is not None:
            self.mu = mu  # Gravitational parameter
        elif m1 is not None:
            self.mu = self._calc_mu(m1, m2)
        else:
            raise ValueError("Provide either mu or m1 (and optionally m2) to calculate mu.")
        
        self.r_vec_0 = r_vec_0
        self.v_vec_0 = v_vec_0
        self.t0 = t0

        self.r0 = np.linalg.norm(r_vec_0)
        self.v0 = np.linalg.norm(v_vec_0)
        self.vr0 = np.dot(r_vec_0, v_vec_0) / self.r0
        self.alpha = 2/self.r0 - self.v0**2/self.mu
    
    def S(self, z):
        """
        Stumpff function S(z).
        """
        if z > 0:
            return (np.sqrt(z) - np.sin(np.sqrt(z)))/np.sqrt(z)**3
        elif z < 0:
            return (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/np.sqrt(-z)**3
        else: # z = 0
            return 1/6

    def C(self, z):
        """
        Stumpff function C(z).
        """
        if z > 0:
            return (1 - np.cos(np.sqrt(z)))/z
        elif z < 0:
            return (np.cosh(np.sqrt(-z)) - 1)/(-z)
        else: # z = 0
            return 1/2

    def state_at_time(self, t):
        """
        Finds the position and velocity vectors at a given time.
        """
        delta_t = t - self.t0
        Q = self.Q_at_delta_t(delta_t)
        r_vec, v_vec = self.state_at_Q(Q, delta_t)
        return r_vec, v_vec

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

    def state_at_Q(self, Q, delta_t):
        """
        Finds the position and velocity vectors at a given universal anomaly and delta time.
        """
        z = self.alpha * Q**2
        
        f = 1 - Q**2/self.r0 * self.C(z)
        g = delta_t - Q**3/np.sqrt(self.mu) * self.S(z)
        
        r_vec = f*self.r_vec_0 + g*self.v_vec_0
        r = np.linalg.norm(r_vec)

        f_dot = -np.sqrt(self.mu)/(r*self.r0) * (self.alpha*Q**3 * self.S(z) - Q)
        g_dot = 1 - Q**2/r * self.C(z)
        v_vec = f_dot*self.r_vec_0 + g_dot*self.v_vec_0

        return r_vec, v_vec
    
