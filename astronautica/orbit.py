from .frame import Frame
from .body import Body
import numpy as np
from scipy.optimize import root_scalar

class Orbit:
    """
    Class representing an orbital system. Provides multiple ways to initialize an orbit:
    - From orbital elements like angular momentum, eccentricity, semi-major axis, inclination,
      right ascension of the ascending node, argument of periapsis and true anomaly
    - From periapsis and apoapsis radii
    - From state vectors (position and velocity)
    ...
    """
    ########### STRING REPRESENTATION ###########
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
            # Format numbers as float with 4 decimal places
            if isinstance(value, (int, float)):
                result += f"{name}: {value:.4f}\n"
            else:
                result += f"{name}: {value}\n"

        return result

    ########### CONSTRUCTORS ###########
    def __init__(self, main_body, _from_classmethod=False):
        """
        Private constructor. Use classmethods to initialize an orbit.
        """
        if not _from_classmethod:
            raise ValueError("Use classmethods like from_elements, from_state_vectors, etc. to initialize an orbit.")
        
        if not isinstance(main_body, Body):
            raise ValueError("main_body must be an instance of Body")
        
        self.points_in_orbital_plane = []
        self.orbital_positions = []
        
        self.main_body = main_body
        self.mu = main_body.mu()
        self.h = 0
        self.e = 0
        self.a = 0
        self.rp = 0
        self.ra = 0
        self.i = 0
        self.omega0 = 0
        self.omega_dot = 0
        self.Omega0 = 0
        self.Omega_dot = 0
        self.theta0 = 0
        self.T = 0
        self.t0_clock = 0
        self.t1_clock = 0
        self.t0_orbit = 0
        self.time_offset = 0
        self.Q_bc_perit0 = np.eye(3)
        self.r_vec_bc = np.zeros(3)
        self.v_vec_bc = np.zeros(3)
        self.h_vec_bc = np.zeros(3)
        self.e_vec_bc = np.zeros(3)
        self.n_vec_bc = np.zeros(3)

    @classmethod
    def from_elements(cls, main_body, h=None, e=None, a=None, rp=None, ra=None, i=0.0, Omega0=0.0, omega0=0.0, theta0=0.0, t0_clock=0.0):
        """
        Initialize an orbit from Keplerian elements.

        :param h: Specific angular momentum (km²/s)
        :param e: Eccentricity
        :param a: Semi-major axis (km)
        :param rp: Periapsis radius (km)
        :param ra: Apoapsis radius (km)
        :param i: Inclination (degrees)
        :param Omega0: Right ascension of the ascending node (degrees)
        :param omega0: Argument of periapsis (degrees)
        :param theta0: True anomaly at t0 (degrees)
        :param t0_clock: Time at t0 (s)
        """
        orbit = cls(main_body, _from_classmethod=True)
        if h is not None and e is not None:
            a = orbit.a_from_h_e(h, e)
            rp = orbit.rp_from_h_e(h, e)
            ra = orbit.ra_from_h_e(h, e)
        elif a is not None and e is not None:
            h = orbit.h_from_a_e(a, e)
            rp = orbit.rp_from_h_e(h, e)
            ra = orbit.ra_from_h_e(h, e)
        elif h is not None and a is not None:
            e = orbit.e_from_h_a(h, a)
            rp = orbit.rp_from_h_e(h, e)
            ra = orbit.ra_from_h_e(h, e)
        elif rp is not None and e is not None:
            h = orbit.h_from_rp_e(rp, e)
            a = orbit.a_from_h_e(h, e)
            ra = orbit.ra_from_h_e(h, e)
        elif rp is not None and ra is not None:
            e = orbit.e_from_rp_ra(rp, ra)
            h = orbit.h_from_rp_e(rp, e)
            a = orbit.a_from_h_e(h, e)
        elif h is not None and rp is not None:
            e = orbit.e_from_h_rp(h, rp)
            a = orbit.a_from_h_e(h, e)
            ra = orbit.ra_from_h_e(h, e)
        else:
            raise ValueError("Invalid combination of parameters. Please provide at least two of h, e, a, rp, ra.")

        orbit.h = h
        orbit.e = e
        orbit.a = a
        orbit.rp = rp
        orbit.ra = ra
        orbit.i = i
        orbit.Omega0 = Omega0
        orbit.omega0 = omega0
        orbit.theta0 = theta0
        orbit.t0_clock = t0_clock
        orbit.finish_init()
        print("Initialized orbit from orbital elements.")
        return orbit

    @classmethod
    def from_apsis(cls, main_body, rp, ra, i=0, Omega0=0, omega0=0, theta0=0, t0_clock=0):
        """
        Initialize an orbit from periapsis and apoapsis radii.

        :param rp: Periapsis radius (km)
        :param ra: Apoapsis radius (km)
        :param i: Inclination (degrees)
        :param Omega0: Right ascension of the ascending node (degrees)
        :param omega0: Argument of periapsis (degrees)
        :param theta0: True anomaly at t0 (degrees)
        :param t0_clock: Time at t0 (s)
        """
        orbit = cls(main_body, _from_classmethod=True)
        e = orbit.e_from_rp_ra(rp, ra)
        h = orbit.h_from_rp_e(rp, e)
        a = orbit.a_from_h_e(h, e)

        orbit.h = h
        orbit.e = e
        orbit.a = a
        orbit.rp = rp
        orbit.ra = ra
        orbit.i = i
        orbit.Omega0 = Omega0
        orbit.omega0 = omega0
        orbit.theta0 = theta0
        orbit.t0_clock = t0_clock
        orbit.finish_init()
        print("Initialized orbit from periapsis and apoapsis radii.")
        return orbit

    @classmethod
    def from_state_vectors(cls, main_body, r_vec_bc, v_vec_bc, t0_clock=0):
        """
        Initialize an orbit from state vectors from a bodycentric reference frame.

        :param r_vec_bc: Position vector in the bodycentric reference frame (km)
        :param v_vec_bc: Velocity vector in the bodycentric reference frame (km/s)
        :param t0_clock: Time at t0 (s)
        """
        r_vec_bc = np.array(r_vec_bc)
        v_vec_bc = np.array(v_vec_bc)
        orbit = cls(main_body, _from_classmethod=True)
        h_vec_bc = cls.h_vec_from_state_vectors(r_vec_bc, v_vec_bc)
        h = np.linalg.norm(h_vec_bc)
        e_vec_bc = orbit.e_vec_from_state_vectors(r_vec_bc, v_vec_bc)
        e = np.linalg.norm(e_vec_bc)
        a = orbit.a_from_h_e(h, e)
        rp = orbit.rp_from_h_e(h, e)
        ra = orbit.ra_from_h_e(h, e)

        i = orbit.i_from_h_vec(h_vec_bc)
        n_vec_bc = cls.n_vec_from_h_vec(h_vec_bc)
        if np.linalg.norm(n_vec_bc) == 0:
            Omega0 = np.degrees(np.arctan2(e_vec_bc[1], e_vec_bc[0]))
        else:
            Omega0 = cls.Omega_from_n_vec(n_vec_bc)
        omega0 = cls.omega_from_e_vec_n_vec(e_vec_bc, n_vec_bc)
        theta0 = cls.theta_from_h_vec_e_vec_r_vec(h_vec_bc, e_vec_bc, r_vec_bc)

        orbit.h = h
        orbit.e = e
        orbit.a = a
        orbit.rp = rp
        orbit.ra = ra
        orbit.i = i
        orbit.Omega0 = Omega0
        orbit.omega0 = omega0
        orbit.theta0 = theta0
        orbit.t0_clock = t0_clock
        orbit.finish_init()
        print("Initialized orbit from state vectors.")
        return orbit

    @classmethod
    def from_r_v_gamma(cls, main_body, r, v, gamma, i=0, Omega0=0, omega0=0, t0_clock=0):
        """
        Creates a new instance of Orbit from distance, velocity and flight path angle.

        :param r: Distance from the primary body center to the position (km)
        :param v: Velocity at the position (km/s)
        :param gamma: Flight path angle (degrees)
        :param i: Inclination (degrees)
        :param Omega0: Right ascension of the ascending node (degrees)
        :param omega0: Argument of periapsis (degrees)
        :param theta0: True anomaly at t0 (degrees)
        :param t0_clock: Time at t0 (s)
        """
        orbit = cls(main_body, _from_classmethod=True)
        gamma_rad = np.radians(gamma)
        v_r = v * np.sin(gamma_rad)
        v_t = v * np.cos(gamma_rad)
        h = r * v_t
        e = np.sqrt(v_r**2*h**2/orbit.mu**2 + (h**2/(orbit.mu*r) - 1)**2)
        a = orbit.a_from_h_e(h, e)
        rp = orbit.rp_from_h_e(h, e)
        ra = orbit.ra_from_h_e(h, e)
        theta0 = orbit.theta_from_r_gamma(r, gamma)
        
        orbit.h = h
        orbit.e = e
        orbit.a = a
        orbit.rp = rp
        orbit.ra = ra
        orbit.i = i
        orbit.Omega0 = Omega0
        orbit.omega0 = omega0
        orbit.theta0 = theta0
        orbit.t0_clock = t0_clock
        orbit.finish_init()
        print("Initialized orbit from distance, velocity and flight path angle.")
        return orbit

    @classmethod
    def from_2_positions(cls, main_body, r1, theta1, r2, theta2, i=0, Omega0=0, omega0=0, t0_clock=0):
        """
        Creates a new instance of Orbit from two orbital positions (r,θ).
        
        :param r1: Radius of the first position (km)
        :param theta1: True anomaly of the first position (degrees)
        :param r2: Radius of the second position (km)
        :param theta2: True anomaly of the second position (degrees)
        :param i: Inclination (degrees)
        :param Omega0: Right ascension of the ascending node (degrees)
        :param omega0: Argument of periapsis (degrees)
        :param theta0: True anomaly at t0 (degrees)
        :param t0_clock: Time at t0 (s)
        """
        orbit = cls(main_body, _from_classmethod=True)
        
        e = orbit.e_from_2_positions(r1, theta1, r2, theta2)
        h = orbit.h_from_e_r_theta(e, r1, theta1)
        a = orbit.a_from_h_e(h, e)
        rp = orbit.rp_from_h_e(h, e)
        ra = orbit.ra_from_h_e(h, e)
        
        orbit.h = h
        orbit.e = e
        orbit.a = a
        orbit.rp = rp
        orbit.ra = ra
        orbit.i = i
        orbit.Omega0 = Omega0
        orbit.omega0 = omega0
        orbit.theta0 = theta2 #using theta2 as the initial position
        orbit.t0_clock = t0_clock
        orbit.finish_init()
        
        # Adding Position 1 (r1, theta1)
        t_clock_1 = orbit.t_clock_at_theta(theta1)
        orbit.add_orbital_position(t_clock_1, name="Position 1")
        
        print("Initialized orbit from two positions.")
        return orbit

    @classmethod
    def init_from_3_vectors(cls, main_body, r1_vec, r2_vec, r3_vec, t0_clock=0):
        """
        Creates a new instance of Orbit from three position vectors using Gibb's method.
        """
        orbit = cls(main_body, _from_classmethod=True)
        r1 = np.linalg.norm(r1_vec)
        r2 = np.linalg.norm(r2_vec)
        r3 = np.linalg.norm(r3_vec)

        C23 = np.cross(r2_vec, r3_vec)

        # Check if the vectors are coplanar
        if np.abs(np.dot(r1_vec/r1, C23/np.linalg.norm(C23))) > 1e-5:
            raise ValueError(f"The vectors are not coplanar. dot product of r1 and C23 is {np.dot(r1_vec/r1, C23/np.linalg.norm(C23))}")
        
        # Calculate N, D, S
        N_vec = r1*(np.cross(r2_vec, r3_vec)) + r2*(np.cross(r3_vec, r1_vec)) + r3*(np.cross(r1_vec, r2_vec))
        D_vec = np.cross(r1_vec, r2_vec) + np.cross(r2_vec, r3_vec) + np.cross(r3_vec, r1_vec)
        S_vec = r1_vec*(r2 - r3) + r2_vec*(r3 - r1) + r3_vec*(r1 - r2)
        N = np.linalg.norm(N_vec)
        D = np.linalg.norm(D_vec)
        
        # Calculate v1, v2, v3
        v1_vec = (np.cross(D_vec, r1_vec) / r2 + S_vec) * np.sqrt(orbit.mu/(N * D))
        v2_vec = (np.cross(D_vec, r2_vec) / r2 + S_vec) * np.sqrt(orbit.mu/(N * D))
        v3_vec = (np.cross(D_vec, r3_vec) / r3 + S_vec) * np.sqrt(orbit.mu/(N * D))

        # Initialize with r2, for more accuracy
        orbit = cls.from_state_vectors(main_body, r2_vec, v2_vec, t0_clock=t0_clock)
        theta_r2 = orbit.theta0
        theta_r1 = orbit.theta_at_state_vectors(r1_vec, v1_vec)
        theta_r3 = orbit.theta_at_state_vectors(r3_vec, v3_vec)
        if theta_r1 > theta_r2:
            theta_r1 = theta_r1 - 360
        if theta_r3 < theta_r2:
            theta_r3 = theta_r3 + 360
        t_clock_r1 = orbit.t_clock_at_theta(theta_r1)
        t_clock_r2 = orbit.t_clock_at_theta(theta_r2)
        t_clock_r3 = orbit.t_clock_at_theta(theta_r3)

        # Replace Position 2 with Position 3, so that the initial position is Position 3
        orbit.theta0 = theta_r3
        r3 = orbit.r_at_theta(theta_r3)
        v3 = orbit.v_at_r(r3)
        orbit.orbital_positions[0]['name'] = "Position 3 - Initial Position"
        orbit.orbital_positions[0]['t_clock'] = t_clock_r3
        orbit.orbital_positions[0]['theta'] = theta_r3
        orbit.orbital_positions[0]['r'] = r3
        orbit.orbital_positions[0]['v'] = v3
        orbit.add_orbital_position(t_clock_r1, name="Position 1")
        orbit.add_orbital_position(t_clock_r2, name="Position 2")
        
        print("Initialized orbit from three position vectors.")
        return orbit

    def finish_init(self):
        self.Omega_dot, self.omega_dot = self.oblateness_correction()

        self.frames = {
            "bodycentric": Frame.bodycentric(),
            "perifocal": Frame(name="perifocal",
                               Omega0_bc_frame=self.Omega0,
                               Omega_dot_bc_frame=self.Omega_dot,
                               omega0_bc_frame=self.omega0,
                               omega_dot_bc_frame=self.omega_dot,
                               i0_bc_frame=self.i,
                               t_clock_bc_frame=self.t0_clock),
            "perifocal_t0": Frame(name="perifocal_t0",
                                  Omega0_bc_frame=self.Omega0,
                                  Omega_dot_bc_frame=0,
                                  omega0_bc_frame=self.omega0,
                                  omega_dot_bc_frame=0,
                                  i0_bc_frame=self.i,
                                  t_clock_bc_frame=self.t0_clock),
            "rotating_bodycentric": Frame.rotating_bodycentric(self.main_body)
        }

        self.T = self.get_T()
        self.t0_orbit = self.t_orbit_at_theta(self.theta0)
        self.time_offset = self.t0_clock - self.t0_orbit
        self.Q_bc_perit0 = self.Q_bc_peri_from_t_clock(self.t0_clock)
        self.r_vec_bc, self.v_vec_bc = self.state_vectors_at_theta(self.theta0, frame="bodycentric")
        self.h_vec_bc = self.h_vec_from_state_vectors(self.r_vec_bc, self.v_vec_bc)
        self.e_vec_bc = self.e_vec_from_state_vectors(self.r_vec_bc, self.v_vec_bc)
        self.n_vec_bc = self.n_vec_from_h_vec(self.h_vec_bc)
        self.h = np.linalg.norm(self.h_vec_bc)
        self.n = np.linalg.norm(self.n_vec_bc)
        self.add_orbital_position(self.t0_clock, name="Initial Position")

    ########### ORBITAL ELEMENTS ###########
    def a_from_h_e(self, h, e):
        """
        Calculates the semi-major axis using specific angular momentum and eccentricity.
        Note: 'a' is negative for hyperbolic orbits and infinite for parabolic orbits.
        """
        if e == 1:
            return float('inf')
        return h**2 / (self.mu * (1 - e**2))

    def e_from_h_a(self, h, a):
        """
        Calculates the eccentricity using specific angular momentum and semi-major axis.
        For hyperbolic orbits (e > 1), 'a' should be negative.
        """
        if a == float('inf'):
            return 1
        return np.sqrt(1 - (h**2 / (self.mu * a)))
    
    def e_from_2_positions(self, r1, theta1, r2, theta2):
        """
        Calculates the eccentricity using two orbital positions.
        """
        theta1_rad = np.radians(theta1)
        theta2_rad = np.radians(theta2)
        return (r2/r1 - 1) / (np.cos(theta1_rad) - (r2/r1)*np.cos(theta2_rad))

    def e_from_h_rp(self, h, rp):
        """
        Calculates the eccentricity using specific angular momentum and periapsis distance.
        """
        return h**2 / (self.mu * rp) - 1

    @staticmethod
    def e_from_rp_ra(rp, ra):
        """
        Calculates the eccentricity using periapsis and apoapsis radii.
        """
        return (ra - rp) / (ra + rp)
    
        """
        Calculates the specific angular momentum using periapsis and eccentricity.
        """
        return np.sqrt(self.mu * rp * (1 + e))

    def e_vec_from_state_vectors(self, r_vec, v_vec):
        """
        Calculates the eccentricity vector using state vectors.
        Both vectors must be in the same reference frame.
        """
        r = np.linalg.norm(r_vec)
        v = np.linalg.norm(v_vec)
        vr = np.dot(r_vec, v_vec)/r
        return 1/self.mu * ( (v**2 - self.mu/r)*r_vec - r*vr*v_vec )

    def h_from_a_e(self, a, e):
        """
        Calculates the specific angular momentum using semi-major axis and eccentricity.
        For hyperbolic orbits (e > 1), 'a' should be negative.
        """
        if e == 1:
            raise ValueError("Can't calculate h from 'a' and 'e' for parabolic (e=1) orbits.")
        return np.sqrt(self.mu * a * (1 - e**2))
    
    def h_from_rp_e(self, rp, e):
        """
        Calculates the specific angular momentum using periapsis and eccentricity.
        """
        return np.sqrt(self.mu * rp * (1 + e))
    
    def h_from_e_r_theta(self, e, r, theta):
        """
        Calculates the specific angular momentum using eccentricity and an orbital position.
        """
        theta_rad = np.radians(theta)
        return np.sqrt(self.mu * r * (1 + e * np.cos(theta_rad)))

    @staticmethod
    def h_vec_from_state_vectors(r_vec, v_vec):
        """
        Calculates the specific angular momentum vector using position and velocity vectors.
        Both vectors must be in the same reference frame.
        """
        return np.cross(r_vec, v_vec)

    def ra_from_h_e(self, h, e):
        """
        Calculates the apoapsis distance using specific angular momentum and eccentricity.
        Note: 'ra' is negative for hyperbolic orbits and infinite for parabolic orbits.
        """
        if e == 1:
            return float('inf')
        else:
            return h**2 / (self.mu * (1 - e))
   
    def rp_from_h_e(self, h, e):
        """
        Calculates the periapsis distance using specific angular momentum and eccentricity.
        """
        return h**2 / (self.mu * (1 + e))
    
    @staticmethod
    def i_from_h_vec(h_vec_bc):
        """
        Calculates the inclination using specific angular momentum vector.
        h vector must be in the bodycentric frame.
        """
        return np.degrees(np.arccos(h_vec_bc[2]/np.linalg.norm(h_vec_bc)))

    @staticmethod
    def n_vec_from_h_vec(h_vec_bc):
        """
        Calculates the node vector using specific angular momentum vector.
        h vector must be in the bodycentric frame.
        """
        return np.cross(np.array([0, 0, 1]), h_vec_bc)
    
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
            theta = np.degrees(np.arctan2(sin_theta, cos_theta))
            if e < 1:
                theta = np.mod(theta, 360)
        return theta


    ########### OTHER ORBITAL PARAMETERS ###########
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

    def get_b(self):
        """
        Calculates the semi-minor axis using semi-major axis and eccentricity.
        Note: For hyperbolic orbits (e > 1), 'a' must be a negative value.
        """
        if self.e < 1:
            return self.a * np.sqrt(1 - self.e**2)
        elif self.e == 1:
            return None
        else: # e > 1
            return - self.a * np.sqrt(self.e**2 - 1)
    
    def get_p(self):
        """
        Calculates the semi-latus rectum using specific angular momentum.
        """
        return self.h**2 / self.mu

    def get_epsilon(self):
        """
        Calculates the specific orbital energy using specific angular momentum and eccentricity.
        """
        return -self.mu**2 / (2 * self.h**2) * (1 - self.e**2)

    def get_C3(self):
        """
        Calculates the hyperbolic energy.
        """
        return self.get_epsilon() * 2

    def get_T(self):
        """
        Calculates the orbital period using gravitational parameter and semi-major axis.
        Note: For hyperbolic orbits (e > 1), 'a' must be a negative value.
        """
        if self.a == float('inf') or self.a < 0:
            return float('inf')
        return 2 * np.pi * np.sqrt(self.a**3 / self.mu)

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
        return self.get_b()

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

    def theta_infinity(self):
        """
        Calculates the true anomaly at infinity for hyperbolic orbits.
        """
        if self.e <= 1:
            raise ValueError("Can't calculate true anomaly at infinity for elliptical or circular orbits.")
        return np.degrees(np.arccos(-1/self.e))
    
    def v_infinity(self):
        """
        Calculates the hyperbolic excess velocity.
        """
        if self.e <= 1:
            raise ValueError("Can't calculate hyperbolic excess velocity for elliptical or parabolic orbits.")
        return np.sqrt(self.mu/self.a)


    ########### DYNAMIC PARAMETERS ###########
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

    def oblateness_correction(self):
        """
        Calculates the oblateness correction factor for node regression and perigee advance.
        """
        if self.e >= 1:
            raise ValueError("Can't calculate oblateness correction for parabolic or hyperbolic orbits.")
        i = np.radians(self.i)
        coeficient = - 3/2 * np.sqrt(self.mu) * self.main_body.J2 * self.main_body.radius**2 / (self.a**(7/2) * (1 - self.e**2)**2)
        Omega_dot = np.degrees(coeficient * np.cos(i))
        omega_dot = np.degrees(coeficient * (5/2 * np.sin(i)**2 - 2)) 
        return Omega_dot, omega_dot

    def r_at_theta(self, theta):
        """
        Calculates the distance from the primary body center to a point on the orbit at a given true anomaly.
        """
        theta = np.radians(theta)
        return self.h**2 / (self.mu) * 1/(1 + self.e * np.cos(theta))
  
    def state_vectors_at_theta(self, theta, frame="bodycentric"):
        """
        Calculates the state vectors at a given true anomaly, in a given frame of reference, at time t0.

        :param theta: True anomaly (degrees)
        :param frame: a Frame object or a string to create a Frame object
        """
        if isinstance(frame, str):
            frame = self.frames[frame]
        
        if frame.name == "bodycentric":
            theta_rad = np.radians(theta)
            r_vec_peri = self.h**2 / (self.mu * (1 + self.e*np.cos(theta_rad))) * np.array([np.cos(theta_rad), np.sin(theta_rad), 0])
            v_vec_peri = self.mu / self.h * np.array([-np.sin(theta_rad), self.e + np.cos(theta_rad), 0])
            if not self.omega_dot and not self.Omega_dot:
                Q_bc_peri = self.Q_bc_perit0
            else:
                t_clock = self.t_orbit_at_theta(theta) + self.time_offset
                Q_bc_peri = self.Q_bc_peri_from_t_clock(t_clock)
            r_vec = Q_bc_peri.T @ r_vec_peri
            v_vec = Q_bc_peri.T @ v_vec_peri
            return r_vec, v_vec
        else:
            t_clock = self.t_orbit_at_theta(theta) + self.time_offset
            r_vec_bc, v_vec_bc = self.state_vectors_at_theta(theta, frame="bodycentric")
            r_vec_frame, v_vec_frame = frame.transform_bc_to_frame(np.column_stack((r_vec_bc, v_vec_bc)), t_clock).T
            return r_vec_frame, v_vec_frame

    def state_vectors_at_t_clock(self, t_clock, frame="bodycentric"):
        """
        Finds the position and velocity vectors at a given time.
        """
        theta = self.theta_at_t_clock(t_clock)
        return self.state_vectors_at_theta(theta, frame=frame)
    
        #######   Time evolution   #######

    def t_clock_at_theta(self, theta):
        """
        Calculates the clock time using the true anomaly.
        """
        t_orbit = self.t_orbit_at_theta(theta)
        return t_orbit + self.time_offset

    def t_orbit_at_theta(self, theta):
        """
        Calculates the orbit time using the true anomaly.
        t_orbit=0 is the time of the perigee passage.
        """
        if self.e < 1: # elliptical
            E = self.E_from_theta(theta)
            M = self.M_from_E(E)
            return self.t_orbit_at_M(M) + self.step_function(theta, 360) * self.T
        elif self.e > 1: # hyperbolic
            E = self.E_from_theta(theta)
            M = self.M_from_E(E)
            return self.t_orbit_at_M(M)
        else: # parabolic
            M = self.M_from_theta(theta)
            return self.t_orbit_at_M(M)

    def theta_at_r(self, r):
        """
        Calculates the true anomaly at a given distance from the primary body center.
        Returns both possible values as a list, with angles between 0 and 360 degrees.
        """
        theta1 = np.arccos((self.h**2 / (self.mu * r) - 1) / self.e)
        theta2 = 2*np.pi - theta1
        return [np.degrees(theta1), np.degrees(theta2)]

    def theta_from_r_gamma(self, r, gamma):
        """
        Calculates the true anomaly using the flight path angle.
        """
        gamma = np.radians(gamma)
        if self.e == 0:  # circular orbit
            theta = 0
        else:
            sin_theta = np.sin(gamma) * (1 + self.h**2/(self.mu*r)) / self.e
            cos_theta = (self.h**2/(self.mu*r) - 1) / self.e
            theta = np.mod(np.degrees(np.arctan2(sin_theta, cos_theta)), 360)
        return theta

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
        e_vec = self.e_vec_from_state_vectors(r_vec, v_vec)
        e = np.linalg.norm(e_vec)
        if e == 0:
            theta = 0
        else:
            cos_theta = np.dot(e_vec, r_vec)/(e * np.linalg.norm(r_vec))
            sin_theta = np.dot(h_vec, np.cross(e_vec, r_vec))/(e * h * np.linalg.norm(r_vec))
            theta = np.mod(np.degrees(np.arctan2(sin_theta, cos_theta)), 360)
        return theta

    def v_escape(self, r):
        """
        Calculates the escape velocity at a given distance.
        """
        return np.sqrt(2 * self.mu / r)

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


    ########### ROTATION MATRICES ###########
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

    def Q_bc_peri_from_t_clock(self, t_clock):
        """
        Calculates the direction cosine matrix from perifocal frame to bodycentric frame at time t0.
        """
        delta_t = t_clock - self.t0_clock
        return self.Q_from_Euler_angles(self.Omega0 + self.Omega_dot * delta_t,
                                        self.i, 
                                        self.omega0 + self.omega_dot * delta_t,
                                        pattern="classical")


    ########### AUXILIARY FUNCTIONS AND VARIABLES ###########
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

        def f_dot(E, M):
            if self.e < 1:
                return 1 - self.e * np.cos(E)
            else:
                return self.e * np.cosh(E) - 1
        
        E = root_scalar(f, fprime=f_dot, x0=E0, args=(M), method='newton').root

        return np.degrees(E)

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
 
    def M_from_theta(self, theta):
        """
        Calculates the mean anomaly using the true anomaly.
        """
        if self.e > 1 or self.e < 1:
            raise ValueError("Can't calculate mean anomaly from true anomaly for elliptical or hyperbolic orbits.")
        
        theta = np.radians(theta)
        return np.degrees(1/6 * np.tan(theta/2)**3 + 1/2 * np.tan(theta/2))

    def E_from_theta(self, theta):
        """
        Calculates the eccentric anomaly using the true anomaly.
        """
        theta = np.radians(theta)
        if self.e < 1:
            return np.degrees(2 * np.arctan(np.sqrt((1 - self.e) / (1 + self.e)) * np.tan(theta/2)))
        elif self.e > 1:
            return np.degrees(2 * np.arctan(np.sqrt((self.e - 1) / (self.e + 1)) * np.tan(theta/2)))
        else:
            raise ValueError("Can't calculate eccentric anomaly from true anomaly for parabolic orbits.")

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


    ########### ADDING POINTS ###########
    def add_orbital_position(self, t_clock=None, theta=None, name=None):
        """
        Adds an orbital position to the orbit at a given clock time.
        """
        if t_clock is None and theta is None:
            raise ValueError("Either t_clock or theta must be provided.")
        if t_clock is None:
            t_clock = self.t_orbit_at_theta(theta) + self.time_offset
        else: # theta is None
            theta = self.theta_at_t_clock(t_clock)
        r = self.r_at_theta(theta)
        v = self.v_at_r(r)
        position = {
            'name': name,
            't_clock': t_clock,
            'theta': theta,
            'r': r,
            'v': v
        }
        
        self.orbital_positions.append(position)
    
    def add_point_in_orbital_plane(self, r, theta, t_clock = 0, name=None):
        """
        Adds a point in the orbital plane at a given distance and true anomaly.
        """
        point = {
            'name': name,
            'r': r,
            'theta': theta,
            't_clock': t_clock
        }
        
        self.points_in_orbital_plane.append(point)

