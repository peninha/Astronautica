from .orbit import Orbit
import numpy as np

class Maneuver:
    """
    Class representing a maneuver.
    """
    def __init__(self, radial, t_or_p, normal, t_clock, name="Maneuver", mode="RTN", orbit_name="Orbit", position_name="Maneuver position"):
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
        self.orbit_name = orbit_name
        self.position_name = position_name

    @classmethod
    def from_RPN(cls, RPN, t_clock, name="Maneuver", orbit_name="Orbit", position_name="Maneuver position"):
        """
        Initialize a maneuver from radial, prograde and normal components.
        """
        return cls(RPN[0], RPN[1], RPN[2], t_clock, name=name, mode="RPN", orbit_name=orbit_name, position_name=position_name)
    
    @classmethod 
    def from_RTN(cls, RTN, t_clock, name="Maneuver", orbit_name="Orbit", position_name="Maneuver position"):
        """
        Initialize a maneuver from radial, tangential and normal components.
        """
        return cls(RTN[0], RTN[1], RTN[2], t_clock, name=name, mode="RTN", orbit_name=orbit_name, position_name=position_name)

    @staticmethod
    def delta_v_for_phase_change(orbit, phase_change, theta_burn=0, n=1, oblateness_correction=True):
        """
        Calculate the delta-v for a phase change in an orbit.
        """
        
        
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

        Parameters
        ----------
        orbit : Orbit
            The orbit to perform the phase change on
        phase_change : float
            The phase change to perform. Positive values are prograde, negative values are retrograde.
        theta_burn : float, optional
            The theta at which the maneuver is performed
        n : int, optional
            The number of orbits to perform the phase change over
        oblateness_correction : bool, optional
            Whether to correct for oblateness
        """
        delta_v1, delta_v2, T1 = cls.delta_v_for_phase_change(orbit, phase_change, theta_burn, n, oblateness_correction)
        phase_maneuver1 = cls.from_RPN([0, delta_v1, 0], orbit.t_clock_at_theta(theta_burn), name="Phase Maneuver 1", orbit_name="Phase Transfer orbit", position_name="Phase starting position")
        phase_maneuver2 = cls.from_RPN([0, delta_v2, 0], orbit.t_clock_at_theta(theta_burn) + T1*n, name="Phase Maneuver 2", orbit_name="Final orbit", position_name="Arrival position")
        return phase_maneuver1, phase_maneuver2
    
    @classmethod
    def orbit_pos_to_orbit_pos_maneuver(cls, orbit_from, orbit_to, theta_burn=0, theta_arrival=0, bounds=None):
        """
        Creates two maneuvers to perform a transfer between two orbits, at specific orbital positions.
        The maneuvers are such that have the least combined delta-v as possible.

        Parameters
        ----------
        orbit_from : Orbit
            The orbit to perform the maneuver from
        orbit_to : Orbit
            The orbit to perform the maneuver to
        theta_burn : float, optional
            The theta at which the first maneuver is performed on orbit_from
        theta_arrival : float, optional
            The theta at which the second maneuver is performed on orbit_to
        bounds : tuple, optional
            The bounds to perform the optimization over
        """
        from scipy.optimize import minimize_scalar

        t_clock_burn = orbit_from.t_clock_at_theta(theta_burn)
        r1_vec, v1_vec = orbit_from.state_vectors_at_theta(theta_burn, frame="bodycentric")
        r2_vec, v2_vec = orbit_to.state_vectors_at_theta(theta_arrival, frame="bodycentric")

        # Objective function to find the optimal delta_t to minimize transfer delta-v
        def objective(delta_t):
            try:
                transfer_orbit = Orbit.from_2_vectors_and_delta_time(
                    orbit_from.main_body, r1_vec, r2_vec, delta_t, t_clock_burn + delta_t)
                
                v_vec_transfer_from = transfer_orbit.state_vectors_at_t_clock(
                    transfer_orbit.orbital_positions[0]['t_clock'], frame="bodycentric")[1]
                delta_v1 = np.linalg.norm(v_vec_transfer_from - v1_vec)
                
                v_vec_transfer_to = transfer_orbit.state_vectors_at_t_clock(
                    transfer_orbit.orbital_positions[1]['t_clock'], frame="bodycentric")[1]
                delta_v2 = np.linalg.norm(v2_vec - v_vec_transfer_to)
                
                delta_v = delta_v1 + delta_v2
                return delta_v
                
            except RuntimeError as e:
                if "failed to converge" in str(e):
                    print(f"Warning: Failed to converge for delta_t={delta_t:.2f}")
                    return np.inf
                raise

        # bounds for searching for delta_t
        if bounds is None:
            bounds = [orbit_from.T/20, orbit_from.T*20]

        # find the optimal delta_t to minimize transfer delta-v
        result = minimize_scalar(
            objective,
            bounds=bounds,
            method='bounded'
        )

        if not result.success:
            raise RuntimeError(f"Failed to converge: {result.message}")

        delta_t = result.x

        t_clock_arrival = t_clock_burn + delta_t

        # Calculate the final delta_v with the optimal time
        transfer_orbit = Orbit.from_2_vectors_and_delta_time(
            orbit_from.main_body, r1_vec, r2_vec, delta_t, t_clock_arrival
        )
        v_vec_transfer_from = transfer_orbit.state_vectors_at_t_clock(t_clock_burn, frame="bodycentric")[1]
        delta_v1_vec = v_vec_transfer_from - v1_vec
        delta_v1 = np.linalg.norm(delta_v1_vec)
        v_vec_transfer_to = transfer_orbit.state_vectors_at_t_clock(t_clock_arrival, frame="bodycentric")[1]
        delta_v2_vec = v2_vec - v_vec_transfer_to
        delta_v2 = np.linalg.norm(delta_v2_vec)
        delta_v = delta_v1 + delta_v2

        # convert the delta_v vectors to RTN components
        RTN1 = transfer_orbit.convert_cartesian_to_RTN(delta_v1_vec, t_clock_burn, frame="bodycentric")
        RTN2 = transfer_orbit.convert_cartesian_to_RTN(delta_v2_vec, t_clock_arrival, frame="bodycentric")

        # create the maneuvers
        maneuver1 = Maneuver.from_RTN(RTN1, t_clock_burn, name="Transfer maneuver", orbit_name="Transfer orbit", position_name="Transfer position")
        maneuver2 = Maneuver.from_RTN(RTN2, t_clock_arrival, name="Orbit insertion maneuver", orbit_name="Final orbit", position_name="Arrival position")

        return maneuver1, maneuver2, delta_t, delta_v
            
    @classmethod
    def hohmann_transfer(cls, orbit_from, orbit_to, theta_burn=0):
        """
        Creates two maneuvers to perform a hohmann transfer between two orbits.
        """
        if np.abs(np.arccos(np.dot(orbit_from.h_vec_bc, orbit_to.h_vec_bc)/(np.linalg.norm(orbit_from.h_vec_bc)*np.linalg.norm(orbit_to.h_vec_bc)))) > 1e-7:
            raise ValueError("Hohmann transfer requires coplanar orbits")

        angle_between_periapsis = np.degrees(np.arctan2(np.linalg.norm(np.cross(orbit_from.e_vec_bc, orbit_to.e_vec_bc)), 
                                           np.dot(orbit_from.e_vec_bc, orbit_to.e_vec_bc)))
        print(angle_between_periapsis)

        theta_arrival = theta_burn + 180 - angle_between_periapsis

        r1_vec, v1_vec = orbit_from.state_vectors_at_theta(theta_burn, frame="bodycentric")
        r2_vec, v2_vec = orbit_to.state_vectors_at_theta(theta_arrival, frame="bodycentric")

        r1 = np.linalg.norm(r1_vec)
        r2 = np.linalg.norm(r2_vec)

        t_clock_burn = orbit_from.t_clock_at_theta(theta_burn)
        transfer_orbit = Orbit.from_2_positions(orbit_from.main_body, r1, 0, r2, 180, i=orbit_from.i, Omega0=orbit_from.Omega0, omega0=theta_burn)

        v_vec_transfer_from = transfer_orbit.state_vectors_at_theta(0, frame="bodycentric")[1]
        delta_v1_vec = v_vec_transfer_from - v1_vec
        delta_v1 = np.linalg.norm(delta_v1_vec)

        v_vec_transfer_to = transfer_orbit.state_vectors_at_theta(180, frame="bodycentric")[1]
        delta_v2_vec = v2_vec - v_vec_transfer_to
        delta_v2 = np.linalg.norm(delta_v2_vec)
        delta_v = delta_v1 + delta_v2

        t_clock_burn = orbit_from.t_clock_at_theta(theta_burn)
        t_clock_arrival = t_clock_burn + transfer_orbit.T/2

        # convert the delta_v vectors to RTN components
        RTN1 = transfer_orbit.convert_cartesian_to_RTN(delta_v1_vec, -transfer_orbit.T/2, frame="bodycentric")
        RTN2 = transfer_orbit.convert_cartesian_to_RTN(delta_v2_vec, 0, frame="bodycentric")

        # create the maneuvers
        maneuver1 = Maneuver.from_RTN(RTN1, t_clock_burn, name="Transfer maneuver", orbit_name="Transfer orbit", position_name="Transfer position")
        maneuver2 = Maneuver.from_RTN(RTN2, t_clock_arrival, name="Orbit insertion maneuver", orbit_name="Final orbit", position_name="Arrival position")

        return maneuver1, maneuver2, delta_v