import numpy as np
import matplotlib.pyplot as plt

class Orbit:
    G = 6.67430e-20  # Constante gravitacional em km^3/kg/s^2
    
    def __init__(self, mu=None, m1=None, m2=0, a=None, e=None, rp=None, ra=None, h=None, epsilon=None, body1radius=None, body2radius=None):
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
    def init_from_points(cls, mu=None, m1=None, m2=0, r1=None, theta1=None, r2=None, theta2=None, body1radius=None, body2radius=None):
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
        # Converter ângulos para radianos
        theta1 = np.radians(theta1)
        theta2 = np.radians(theta2)
        
        # Calcular e
        e = (r2/r1 - 1) / (np.cos(theta1) - (r2/r1)*np.cos(theta2))
        
        # Calcular h
        if mu is None and m1 is not None:
            mu = cls.G * (m1 + m2)
            
        if mu is None:
            raise ValueError("Forneça mu ou m1 (e opcionalmente m2) para calcular mu.")
            
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
        gamma = np.radians(gamma)
        vtan = v * np.sin(gamma)
        vrad = v * np.cos(gamma)
        h = r * vrad
        e = vtan**2 / mu - 1/r
        return cls(mu=mu, e=e, h=h, body1radius=body1radius, body2radius=body2radius)

    ############# Calculations #############

    #####   mu   #####
    def _calc_mu(self, m1, m2=0):
        """
        Calculates the gravitational parameter mu based on the masses of two bodies.
        
        :param m1: Mass of the primary body (kg)
        :param m2: Mass of the secondary body, default is 0 (kg)
        """
        return self.G * (m1 + m2)
  
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

    def r_mean(self):
        """
        Calculates the mean radius using periapsis and apoapsis distances.
        Note: For hyperbolic orbits (e > 1), 'ra' must be a negative value.
        """
        return (self._rp * self._ra)**(1/2)

    def v_at_r(self, r):
        """
        Calculates the orbital velocity at a given distance using the vis-viva equation.
        
        :param r: Distance from the primary body center to the point (km)
        :return: Velocity at the point (km/s)
        """
        return np.sqrt(self.mu * (2 / r - 1 / self._a))
    
    def r_at_theta(self, theta):
        """
        Calculates the distance from the primary body center to a point on the orbit at a given true anomaly.
        """
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

    def plot(self, points=None):
        """
        Plots the orbit with optional points.
        """
        # Criar array de anomalia verdadeira
        if self._e == 1:
            theta = np.linspace(-np.pi*2/3, np.pi*2/3, 1000)
        elif self._e > 1:
            theta = np.linspace(-self.theta_inf(), self.theta_inf(), 1000)
        else:
            theta = np.linspace(0, 2*np.pi, 1000)
        
        # Calcular raio para cada ângulo
        r = self.r_at_theta(theta)
        
        # Converter para coordenadas cartesianas
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        
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

        
