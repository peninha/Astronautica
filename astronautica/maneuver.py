# astronautica/maneuver.py
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
    def delta_v_for_phase_change(orbit, phase_change, theta_burn=0, n=1):
        """
        Calculate the delta-v for a phase change in an orbit.
        """
        v0 = orbit.v_at_theta(theta_burn)
        a0 = orbit.a
        theta_prime = theta_burn + phase_change
        T1 = orbit.t_clock_at_theta(theta_prime + 360 - phase_change) - orbit.t_clock_at_theta(theta_prime)
        a1 = (T1 * np.sqrt(orbit.mu)/(2*np.pi))**(2/3)
        rp1 = orbit.rp
        ra1 = 2*a1 - rp1
        e1 = orbit.e_from_rp_ra(rp1, ra1)
        h1 = orbit.h_from_a_e(a1, e1)
        v1 = h1 / rp1
        delta_v = v1 - v0
        return delta_v, -delta_v

    