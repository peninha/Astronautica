import numpy as np

class Maneuver:
    """
    Class representing a maneuver.
    """
    def __init__(self, radial, t_or_p, normal, t_clock, name=None, mode="RTN"):
        """
        Initialize a maneuver.

        Parameters
        ----------
        radial : float
            Radial component of the maneuver
        prograde : float 
            Prograde/tangential component of the maneuver
        normal : float
            Normal component of the maneuver
        t_clock : float
            Time of the maneuver
        name : str, optional
            Name of the maneuver
        mode : str, optional
            Mode of initialization - "RTN" or "RPN"
        """
        if mode == "RTN":
            self.RTN = np.array([radial, t_or_p, normal])
            self.RPN = None
        elif mode == "RPN":
            self.RPN = np.array([radial, t_or_p, normal]) 
            self.RTN = None
        else:
            raise ValueError("Mode must be 'RTN' or 'RPN'")
        
        self.mode = mode
        self.t_clock = t_clock
        self.name = name
        self.delta_v_vec_bc = None

    @classmethod
    def from_RPN(cls, radial, prograde, normal, t_clock, name=None):
        """
        Initialize a maneuver from radial, prograde and normal components.
        """
        return cls(radial, prograde, normal, t_clock, name, mode="RPN")

    @classmethod 
    def from_RTN(cls, radial, tangential, normal, t_clock, name=None):
        """
        Initialize a maneuver from radial, tangential and normal components.
        """
        return cls(radial, tangential, normal, t_clock, name, mode="RTN")

    @staticmethod
    def delta_v_for_phase_change(orbit, phase_change, theta_burn=0, n=1, oblateness_correction=True):
        """
        Calculate the delta-v for a phase change in an orbit.
        """
        from astronautica.orbit import Orbit
        
        v0 = orbit.v_at_theta(theta_burn)
        r = orbit.r_at_theta(theta_burn)
        t_clock_burn = orbit.t_clock_at_theta(theta_burn)
        orbit_target = Orbit.from_elements(orbit.main_body, h=orbit.h, e=orbit.e, theta0=theta_burn + phase_change, i=orbit.i, Omega0=orbit.Omega0, omega0=orbit.omega0, t0_clock=t_clock_burn)
        
        T_burn_to_phase = orbit.t_clock_at_theta(theta_burn + phase_change) - orbit.t_clock_at_theta(theta_burn)
        delta_T = T_burn_to_phase/n
        T1 = orbit.T - delta_T

        ### Correcting for oblateness ###
        def calculate_delta_theta():
            v1 = np.sqrt((2*orbit.mu/r) - (2*np.pi*orbit.mu/T1)**(2/3))
            r_vec, v_vec = orbit.state_vectors_at_t_clock(t_clock_burn, frame="bodycentric")
            v = np.linalg.norm(v_vec)
            v1_vec = v_vec*v1/v
            orbit1 = Orbit.from_state_vectors(orbit.main_body, r_vec, v1_vec, t0_clock=t_clock_burn)
            r_vec1 = orbit1.state_vectors_at_t_clock(t_clock_burn + T1*n, frame="bodycentric")[0]
            r_vec_target = orbit_target.state_vectors_at_t_clock(t_clock_burn + T1*n, frame="bodycentric")[0]
            theta_end_1 = orbit1.convert_cartesian_to_polar(r_vec1)[1]
            theta_end_target = orbit_target.convert_cartesian_to_polar(r_vec_target)[1]
            return theta_end_1 - theta_end_target
        
        if oblateness_correction:
            delta_theta = calculate_delta_theta()
            while np.abs(delta_theta) > 0.5:
                extra_time = (orbit.t_clock_at_theta(theta_burn + delta_theta) - orbit.t_clock_at_theta(theta_burn))/n
                T1 += extra_time
                delta_theta = calculate_delta_theta()
        
        v1 = np.sqrt((2*orbit.mu/r) - (2*np.pi*orbit.mu/T1)**(2/3))      
        delta_v = v1 - v0
        return delta_v, -delta_v, T1

    @classmethod
    def phase_maneuver(cls, orbit, phase_change, theta_burn=0, n=1, oblateness_correction=True):
        """
        Calculate the delta-v for a phase change in an orbit.
        """
        delta_v1, delta_v2, T1 = cls.delta_v_for_phase_change(orbit, phase_change, theta_burn, n, oblateness_correction)
        phase_maneuver1 = cls.from_RPN(0, delta_v1, 0, orbit.t_clock_at_theta(theta_burn), name="Phase Maneuver 1")
        phase_maneuver2 = cls.from_RPN(0, delta_v2, 0, orbit.t_clock_at_theta(theta_burn) + T1*n, name="Phase Maneuver 2")
        return phase_maneuver1, phase_maneuver2