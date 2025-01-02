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
        T_burn_to_phase = orbit.t_clock_at_theta(theta_burn + phase_change) - orbit.t_clock_at_theta(theta_burn)
        delta_T = T_burn_to_phase/n
        T1 = orbit.T - delta_T
        v0 = orbit.v_at_theta(theta_burn)
        r = orbit.r_at_theta(theta_burn)
        v1 = np.sqrt((2*orbit.mu/r) - (2*np.pi*orbit.mu/T1)**(2/3))
        delta_v = v1 - v0
        return delta_v, -delta_v, T1

    @classmethod
    def phase_maneuver(cls, orbit, phase_change, theta_burn=0, n=1):
        """
        Calculate the delta-v for a phase change in an orbit.
        """
        delta_v1, delta_v2, T1 = cls.delta_v_for_phase_change(orbit, phase_change, theta_burn, n)
        phase_maneuver1 = cls.from_RPN(0, delta_v1, 0, orbit.t_clock_at_theta(theta_burn), name="Phase Maneuver 1")
        phase_maneuver2 = cls.from_RPN(0, delta_v2, 0, orbit.t_clock_at_theta(theta_burn) + T1*n, name="Phase Maneuver 2")
        return phase_maneuver1, phase_maneuver2