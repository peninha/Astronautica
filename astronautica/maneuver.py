from aiohttp import ClientSession
from .orbit import Orbit
import numpy as np
from scipy.optimize import minimize_scalar, minimize, root_scalar
import math

class Maneuver:
    """
    Class representing a maneuver.
    """
    def __init__(self, RXN, t_clock, name="Maneuver", mode="RTN", orbit_name="Orbit", position_name="Maneuver position"):
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
            self.RTN = np.array(RXN)
            self.RPN = None
            self.delta_v = np.linalg.norm(self.RTN)
        elif mode == "RPN":
            self.RPN = np.array(RXN) 
            self.RTN = None
            self.delta_v = None
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
        return cls(RPN, t_clock, name=name, mode="RPN", orbit_name=orbit_name, position_name=position_name)
    
    @classmethod 
    def from_RTN(cls, RTN, t_clock, name="Maneuver", orbit_name="Orbit", position_name="Maneuver position"):
        """
        Initialize a maneuver from radial, tangential and normal components.
        """
        return cls(RTN, t_clock, name=name, mode="RTN", orbit_name=orbit_name, position_name=position_name)

    @classmethod
    def from_delta_v_and_angle(cls, delta_v, angle, t_clock, name="Maneuver", orbit_name="Orbit", position_name="Maneuver position"):
        """
        Initialize a maneuver from delta-v and angle to local horizon.

        Parameters
        ----------
        delta_v : float
            Delta-v of the maneuver
        angle : float
            Angle to local horizon
        """
        angle = np.radians(angle)
        RPN = np.array([delta_v*np.sin(angle), delta_v*np.cos(angle), 0])
        return cls(RPN, t_clock, name=name, mode="RPN", orbit_name=orbit_name, position_name=position_name)

    @classmethod
    def from_delta_v_vec(cls, delta_v_vec, orbit, t_clock, name="Maneuver", orbit_name="Orbit", position_name="Maneuver position"):
        """
        Initialize a maneuver from delta-v vector and theta.
        """
        RTN = orbit.convert_cartesian_to_RTN(delta_v_vec, t_clock, frame="bodycentric")
        return cls(RTN, t_clock, name=name, mode="RTN", orbit_name=orbit_name, position_name=position_name)
    
    @staticmethod
    def delta_v_for_phase_change(orbit, phase_change, theta_burn=0.0, n=1, oblateness_correction=True):
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
    def phase_maneuver(cls, orbit, phase_change, theta_burn=0.0, n=1, oblateness_correction=True):
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
    
    @staticmethod
    def orbit_pos_to_orbit_pos_maneuver(orbit_from, orbit_to, theta_burn=0.0, theta_arrival=0.0, bounds=None, max_iter=1000, tol=1e-6):
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
            method='bounded',
            options={'xatol': tol, 'maxiter': max_iter}
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
            
    @staticmethod
    def hohmann_transfer(orbit_from, orbit_to, theta_burn=0.0):
        """
        Creates two maneuvers to perform a hohmann transfer between two orbits.
        """
        if np.abs(np.arccos(np.dot(orbit_from.h_vec_bc, orbit_to.h_vec_bc)/(np.linalg.norm(orbit_from.h_vec_bc)*np.linalg.norm(orbit_to.h_vec_bc)))) > 1e-7:
            raise ValueError("Hohmann transfer requires coplanar orbits")

        if orbit_to.e < 1e-7:
            theta_arrival = theta_burn + 180
            r1_vec, v1_vec = orbit_from.state_vectors_at_theta(theta_burn, frame="bodycentric")
            r1 = np.linalg.norm(r1_vec)
            r2 = orbit_to.r_at_theta(0)
            r2_vec = orbit_from.state_vectors_at_theta(theta_arrival, frame="bodycentric")[0]/r1 * r2
            v2_vec_peri = np.sqrt(orbit_to.mu/r2) * np.array([-np.sin(np.radians(theta_arrival)), np.cos(np.radians(theta_arrival)), 0])
            v2_vec = orbit_from.Q_bc_peri_from_t_clock(0).T @ v2_vec_peri

        else:
            angle_between_periapsis = np.degrees(np.arctan2(np.linalg.norm(np.cross(orbit_from.e_vec_bc, orbit_to.e_vec_bc)), 
                                           np.dot(orbit_from.e_vec_bc, orbit_to.e_vec_bc)))
            theta_arrival = theta_burn + 180 - angle_between_periapsis
            r1_vec, v1_vec = orbit_from.state_vectors_at_theta(theta_burn, frame="bodycentric")
            r2_vec, v2_vec = orbit_to.state_vectors_at_theta(theta_arrival, frame="bodycentric")
            r1 = np.linalg.norm(r1_vec)
            r2 = np.linalg.norm(r2_vec)

        t_clock_burn = orbit_from.t_clock_at_theta(theta_burn)
        transfer_orbit = Orbit.from_2_positions(orbit_from.main_body, r1, 0, r2, 180, i=orbit_from.i, Omega0=orbit_from.Omega0, omega0=theta_burn+orbit_from.omega0)

        v_vec_transfer_from = transfer_orbit.state_vectors_at_theta(0, frame="bodycentric")[1]
        delta_v1_vec = v_vec_transfer_from - v1_vec
        delta_v1 = np.linalg.norm(delta_v1_vec)

        v_vec_transfer_to = transfer_orbit.state_vectors_at_theta(180, frame="bodycentric")[1]
        delta_v2_vec = v2_vec - v_vec_transfer_to
        delta_v2 = np.linalg.norm(delta_v2_vec)
        delta_v = delta_v1 + delta_v2

        t_clock_arrival = t_clock_burn + transfer_orbit.T/2

        # convert the delta_v vectors to RTN components
        RTN1 = transfer_orbit.convert_cartesian_to_RTN(delta_v1_vec, -transfer_orbit.T/2, frame="bodycentric")
        RTN2 = transfer_orbit.convert_cartesian_to_RTN(delta_v2_vec, 0, frame="bodycentric")

        # create the maneuvers
        maneuver1 = Maneuver.from_RTN(RTN1, t_clock_burn, name="Transfer maneuver", orbit_name="Transfer orbit", position_name="Transfer position")
        maneuver2 = Maneuver.from_RTN(RTN2, t_clock_arrival, name="Orbit insertion maneuver", orbit_name="Final orbit", position_name="Arrival position")

        return maneuver1, maneuver2, delta_v
    
    @staticmethod
    def bielliptic_transfer(orbit_from, orbit_to, theta_burn=0.0, apoapsis=None):
        """
        Creates two maneuvers to perform a bielliptic transfer between two orbits.

        Parameters
        ----------
        apoapsis : float, optional
            The apoapsis of the transfer orbit
        """
        if apoapsis is None:
            raise ValueError("Apoapsis must be specified for bielliptic transfer")
        
        if np.abs(np.arccos(np.dot(orbit_from.h_vec_bc, orbit_to.h_vec_bc)/(np.linalg.norm(orbit_from.h_vec_bc)*np.linalg.norm(orbit_to.h_vec_bc)))) > 1e-7:
            raise ValueError("Bielliptic transfer requires coplanar orbits")

        if orbit_to.e < 1e-7:
            theta_arrival = theta_burn
            r1_vec, v1_vec = orbit_from.state_vectors_at_theta(theta_burn, frame="bodycentric")
            r1 = np.linalg.norm(r1_vec)
            r3 = orbit_to.r_at_theta(0)
            r3_vec = orbit_from.state_vectors_at_theta(theta_arrival, frame="bodycentric")[0]/r1 * r3
            v3_vec_peri = np.sqrt(orbit_to.mu/r3) * np.array([-np.sin(np.radians(theta_arrival)), np.cos(np.radians(theta_arrival)), 0])
            v3_vec = orbit_from.Q_bc_peri_from_t_clock(0).T @ v3_vec_peri

        else:
            angle_between_periapsis = np.degrees(np.arctan2(np.linalg.norm(np.cross(orbit_from.e_vec_bc, orbit_to.e_vec_bc)), 
                                           np.dot(orbit_from.e_vec_bc, orbit_to.e_vec_bc)))
            theta_arrival = theta_burn - angle_between_periapsis
            r1_vec, v1_vec = orbit_from.state_vectors_at_theta(theta_burn, frame="bodycentric")
            r3_vec, v3_vec = orbit_to.state_vectors_at_theta(theta_arrival, frame="bodycentric")
            r1 = np.linalg.norm(r1_vec)
            r3 = np.linalg.norm(r3_vec)

        r2 = apoapsis
        if r3/r1 < 11.94:
            print("Warning: r3/r1 =", r3/r1, ".Bielliptic transfer is less efficient than hohmann transfer.")
        
        t_clock_burn = orbit_from.t_clock_at_theta(theta_burn)

        transfer_orbit1 = Orbit.from_2_positions(orbit_from.main_body, r1, 0, r2, 180, i=orbit_from.i, Omega0=orbit_from.Omega0, omega0=theta_burn+orbit_from.omega0)
        v_vec_transfer1_from = transfer_orbit1.state_vectors_at_theta(0, frame="bodycentric")[1]
        delta_v1_vec = v_vec_transfer1_from - v1_vec
        delta_v1 = np.linalg.norm(delta_v1_vec)

        transfer_orbit2 = Orbit.from_2_positions(orbit_from.main_body, r2, 180, r3, 0, i=orbit_from.i, Omega0=orbit_from.Omega0, omega0=theta_burn+orbit_from.omega0)
        v_vec_transfer1_apoapsis = transfer_orbit1.state_vectors_at_theta(180, frame="bodycentric")[1]
        v_vec_transfer2_apoapsis = transfer_orbit2.state_vectors_at_theta(180, frame="bodycentric")[1]
        delta_v2_vec = v_vec_transfer2_apoapsis - v_vec_transfer1_apoapsis
        delta_v2 = np.linalg.norm(delta_v2_vec)
       
        v_vec_transfer2_to = transfer_orbit2.state_vectors_at_theta(0, frame="bodycentric")[1]
        delta_v3_vec = v3_vec - v_vec_transfer2_to
        delta_v3 = np.linalg.norm(delta_v3_vec)
        
        delta_v = delta_v1 + delta_v2 + delta_v3

        t_clock_transfer_apoapsis = t_clock_burn + transfer_orbit1.T/2
        t_clock_arrival = t_clock_transfer_apoapsis + transfer_orbit2.T/2

        # convert the delta_v vectors to RTN components
        RTN1 = transfer_orbit1.convert_cartesian_to_RTN(delta_v1_vec, -transfer_orbit1.T/2, frame="bodycentric")
        RTN2 = transfer_orbit2.convert_cartesian_to_RTN(delta_v2_vec, -transfer_orbit2.T/2, frame="bodycentric")
        RTN3 = transfer_orbit2.convert_cartesian_to_RTN(delta_v3_vec, 0, frame="bodycentric")

        # create the maneuvers
        maneuver1 = Maneuver.from_RTN(RTN1, t_clock_burn, name="Transfer maneuver1", orbit_name="Transfer orbit1", position_name="Transfer position")
        maneuver2 = Maneuver.from_RTN(RTN2, t_clock_transfer_apoapsis, name="Transfer maneuver2", orbit_name="Transfer orbit2", position_name="Transfer position")
        maneuver3 = Maneuver.from_RTN(RTN3, t_clock_arrival, name="Orbit insertion maneuver", orbit_name="Final orbit", position_name="Arrival position")

        return maneuver1, maneuver2, maneuver3, delta_v
    
    @staticmethod
    def check_orbit_intersection(orbit_from, orbit_to, theta_burn=0.0, tol=1e-3, num_initial_guesses=12):
        """
        Check if there is an intersection between two orbits at a given theta_burn.

        Parameters
        ----------
        orbit_from : Orbit
            Initial orbit
        orbit_to : Orbit
            Final orbit
        theta_burn : float
            Angle theta at the initial orbit where the intersection is checked
        tol : float, optional
            Tolerance to consider two points as a intersection (in km)
        num_initial_guesses : int, optional
            Number of initial points to search for intersections

        Returns
        -------
        tuple
            (bool, float, np.array) - (exists_intersection, intersection_theta, delta_v_vec)
            If there is no intersection, returns (False, None, None)
        """

        # Get the position vector at theta_burn
        r_burn_vec = orbit_from.state_vectors_at_theta(theta_burn, frame="bodycentric")[0]

        def distance_function(theta2):
            r2_vec = orbit_to.state_vectors_at_theta(theta2[0], frame="bodycentric")[0]
            return np.linalg.norm(r_burn_vec - r2_vec)
        
        # Try different initial points uniformly distributed
        thetas2 = np.linspace(0, 360, num_initial_guesses)
        min_distance = float('inf')
        best_theta2 = None
        
        for theta2 in thetas2:
            result = minimize(
                distance_function,
                x0=theta2,
                method='Nelder-Mead',
                options={'xatol': 1e-7}
            )
            
            if result.fun < min_distance:
                min_distance = result.fun
                best_theta2 = result.x[0] % 360

        # If an intersection is found within the tolerance
        if min_distance < tol:
            # Calculate the delta_v required
            v1_vec = orbit_from.state_vectors_at_theta(theta_burn, frame="bodycentric")[1]
            v2_vec = orbit_to.state_vectors_at_theta(best_theta2, frame="bodycentric")[1]
            delta_v_vec = v2_vec - v1_vec
            return True, best_theta2, delta_v_vec
        
        return False, None, None

    @staticmethod
    def find_orbit_intersections(orbit1, orbit2, tol=1e-3, num_initial_guesses=12):
        """
        Find the points where two orbits intersect using numerical optimization.

        Parameters
        ----------
        orbit1 : Orbit
            First orbit
        orbit2 : Orbit
            Second orbit
        tol : float, optional
            Tolerance to consider two points as a intersection (in km)
        num_initial_guesses : int, optional
            Number of initial points to search for intersections

        Returns
        -------
        list
            List of dictionaries containing thetas and positions of intersection points
        """
        def distance_function(x):
            theta1, theta2 = x
            r1 = orbit1.state_vectors_at_theta(theta1, frame="bodycentric")[0]
            r2 = orbit2.state_vectors_at_theta(theta2, frame="bodycentric")[0]
            return np.linalg.norm(r1 - r2)
        
        intersections = []
        # Generate pairs of initial angles uniformly distributed
        thetas1 = np.linspace(0, 360, num_initial_guesses)
        thetas2 = np.linspace(0, 360, num_initial_guesses)
        
        for theta1 in thetas1:
            for theta2 in thetas2:
                # Optimization to find exact intersection point
                result = minimize(
                    distance_function,
                    x0=[theta1, theta2],
                    method='Nelder-Mead',
                    options={'xatol': 1e-7}
                )
                
                if result.fun < tol:
                    theta1_opt, theta2_opt = result.x
                    position = orbit1.state_vectors_at_theta(theta1_opt, frame="bodycentric")[0]
                    
                    # Check if this point was already found
                    is_new = True
                    for known in intersections:
                        if np.linalg.norm(known['position'] - position) < tol:
                            is_new = False
                            break
                    
                    if is_new:
                        intersections.append({
                            'theta1': theta1_opt % 360,
                            'theta2': theta2_opt % 360,
                            'position': position
                        })
        
        return intersections

    @classmethod
    def orbit_change_maneuver_at_theta(cls, orbit_from, orbit_to, theta_burn, tol=1e-3, num_initial_guesses=12):
        """
        Creates a maneuver to change the orbit of a spacecraft at a given theta_burn relative to orbit_from.
        This theta_burn must be an intersection point between the two orbits.

        Parameters
        ----------
        orbit_from : Orbit
            Initial orbit
        orbit_to : Orbit
            Final orbit
        theta_burn : float
            The theta_burn at which the maneuver is performed
        tol : float, optional
            Tolerance to consider two points as a intersection (in km)
        num_initial_guesses : int, optional
            Number of initial points to search for intersections

        Returns
        -------
        tuple
            (maneuver, theta_burn, delta_v)
        """
        
        intersect, _ , delta_v_vec = cls.check_orbit_intersection(orbit_from, orbit_to, theta_burn, tol, num_initial_guesses)
        if not intersect:
            raise ValueError("The orbits do not intersect within the specified tolerance")
        
        t_clock_burn = orbit_from.t_clock_at_theta(theta_burn)
        RTN = orbit_from.convert_cartesian_to_RTN(delta_v_vec, t_clock_burn, frame="bodycentric")
        maneuver = cls.from_RTN(RTN, t_clock_burn, name="Change orbit maneuver", orbit_name="Final orbit", position_name="Intersection point")
        delta_v = float(np.linalg.norm(delta_v_vec)) # type: ignore
        return maneuver, theta_burn, delta_v

    @classmethod
    def orbit_change_maneuver(cls, orbit_from, orbit_to, tol=1e-3, num_initial_guesses=12):
        """
        Creates a maneuver to change the orbit of a spacecraft.
        The two orbits must intersect at least at one point.
        The maneuver is performed at the intersection point that requires the least delta-v.

        Parameters
        ----------
        orbit_from : Orbit
            Initial orbit
        orbit_to : Orbit
            Final orbit

        Returns
        -------
        tuple
            (maneuver, theta_burn, delta_v)
        """
        intersections = cls.find_orbit_intersections(orbit_from, orbit_to, tol, num_initial_guesses)
        
        if not intersections:
            raise ValueError("The orbits do not intersect within the specified tolerance")
        
        # Find the intersection point with the least delta-v
        min_delta_v = float('inf')
        best_maneuver = None
        chosen_theta_burn = None
        
        for intersection in intersections:
            theta_burn = intersection['theta1']
            theta_arrival = intersection['theta2']
            
            v1_vec = orbit_from.state_vectors_at_theta(theta_burn, frame="bodycentric")[1]
            v2_vec = orbit_to.state_vectors_at_theta(theta_arrival, frame="bodycentric")[1]
            
            delta_v_vec = v2_vec - v1_vec
            delta_v = np.linalg.norm(delta_v_vec)
            
            if delta_v < min_delta_v:
                min_delta_v = delta_v
                chosen_theta_burn = theta_burn
                t_clock_burn = orbit_from.t_clock_at_theta(theta_burn)
                RTN = orbit_from.convert_cartesian_to_RTN(delta_v_vec, t_clock_burn, frame="bodycentric")
                best_maneuver = cls.from_RTN(
                    RTN, 
                    t_clock_burn, 
                    name="Change orbit maneuver",
                    orbit_name="Final orbit",
                    position_name="Intersection point"
                )
        return best_maneuver, chosen_theta_burn, min_delta_v
    
    @classmethod
    def impact_maneuver_from_t_clock_burn(cls, orbit_from, orbit_to, delta_v_target, delta_t_guess, delta_t_min, delta_t_max, t_clock_burn=0, tol=1e-6, max_iter=1000):
        """
        Finds a maneuver, consuming the whole given delta-v, to impact the target.
        This functions iterates over delta_t to find the maneuver which spend the given delta_v_target.

        Parameters
        ----------
        delta_v_target : float
            Target delta-v for the maneuver
        delta_t_guess : float
            Initial guess for the delta-t
        delta_t_min : float
            Minimum delta-t for the maneuver
        delta_t_max : float
            Maximum delta-t for the maneuver
        t_clock_burn : float, optional
            t_clock at which the maneuver is performed
        tol : float, optional
            Tolerance for the optimization
        max_iter : int, optional
            Maximum number of iterations for the optimization
        """
        def objective(delta_t):
            if isinstance(delta_t, (list, np.ndarray)):
                delta_t = delta_t[0]
            t_clock_impact = t_clock_burn + delta_t
            r_vec_from = orbit_from.state_vectors_at_t_clock(t_clock_burn)[0]
            r_vec_impact = orbit_to.state_vectors_at_t_clock(t_clock_impact)[0]
            impact_orbit = Orbit.from_2_vectors_and_delta_time(orbit_from.main_body, r_vec_from, r_vec_impact, delta_t, t0_clock=t_clock_impact, max_iter=max_iter, tol=tol)
            delta_v = cls.delta_v_from_v_vecs(orbit_from.state_vectors_at_t_clock(t_clock_burn)[1], impact_orbit.state_vectors_at_t_clock(t_clock_burn)[1])
            return delta_v - delta_v_target

        result = root_scalar(objective,
                            bracket=[delta_t_min, delta_t_max],
                            x0=delta_t_guess,
                            method='newton',
                            options={'xtol': tol, 'maxiter': max_iter})

        if result.converged:
            delta_t_opt = result.root
        else:
            raise RuntimeError("Failed to find a solution")

        t_clock_impact = t_clock_burn + delta_t_opt
        r_vec_from = orbit_from.state_vectors_at_t_clock(t_clock_burn)[0]
        r_vec_impact = orbit_to.state_vectors_at_t_clock(t_clock_impact)[0]
        impact_orbit = Orbit.from_2_vectors_and_delta_time(orbit_from.main_body, r_vec_from, r_vec_impact, delta_t_opt, t0_clock=t_clock_impact)

        theta_burn = orbit_from.theta_at_t_clock(t_clock_burn)
        maneuver, _ , delta_v = cls.orbit_change_maneuver_at_theta(orbit_from, impact_orbit, theta_burn)
        return maneuver, delta_t_opt, delta_v

    @classmethod
    def impact_maneuver(cls, orbit_from, orbit_to, theta_rel, delta_v_target, delta_t_guess, delta_t_min, delta_t_max, from_t_clock=300, tol=1e-6, max_iter=1000):
        """
        Finds a maneuver using a given delta-v to impact the target.

        Parameters
        ----------
        theta_rel : float
            Relative angle between the two orbits. Positive if orbit_to is ahead, negative if behind.
        delta_v_target : float
            Target delta-v for the maneuver
        delta_t_min : float
            Minimum delta-t for the maneuver
        delta_t_max : float
            Maximum delta-t for the maneuver
        from_t_clock : float, optional
            t_clock from which the theta_rel is searched, to allow time for maneuver preparation. Default is 300 seconds (5 minutes).
        tolerance : float, optional
            Tolerance for the optimization
        max_iterations : int, optional
            Maximum number of iterations for the optimization
        """
        t_clock_burn = cls.find_t_clock_for_relative_theta(orbit_from, orbit_to, theta_rel, from_t_clock=from_t_clock)
        maneuver, delta_t_opt, delta_v = cls.impact_maneuver_from_t_clock_burn(orbit_from, orbit_to, delta_v_target, delta_t_guess, delta_t_min, delta_t_max, t_clock_burn=t_clock_burn, tol=tol, max_iter=max_iter)
        return {
            "maneuver": maneuver,
            "t_clock_burn": t_clock_burn,
            "delta_t": delta_t_opt,
            "delta_v": delta_v
        }

    @staticmethod
    def relative_periapsis_argument(orbit_from, orbit_to):
        """
        Finds the relative periapsis argument between two coplanar (or nearly coplanar) orbits.
        """
        r_vec_from = orbit_from.state_vectors_at_theta(0, frame="bodycentric")[0]
        r_vec_to = orbit_to.state_vectors_at_theta(0, frame="bodycentric")[0]
        r_from = np.linalg.norm(r_vec_from)
        r_to = np.linalg.norm(r_vec_to)
        omega_rel = np.degrees(np.arccos(np.dot(r_vec_from, r_vec_to)/(r_from*r_to)))
        if np.dot(np.cross(r_vec_from, r_vec_to), orbit_from.h_vec_bc) < 0:
            omega_rel *= -1
        return omega_rel
    
    @classmethod
    def relative_theta0(cls, orbit_from, orbit_to):
        """
        Finds the relative theta0 between two coplanar (or nearly coplanar) orbits.
        """
        omega_rel = cls.relative_periapsis_argument(orbit_from, orbit_to)
        theta0_rel = orbit_to.theta0 - orbit_from.theta0 + omega_rel
        return theta0_rel
    
    @classmethod
    def relative_theta_at_t_clock(cls, orbit_from, orbit_to, t_clock):
        """
        Finds the relative theta at a given time between two coplanar (or nearly coplanar) orbits.
        """
        omega_rel = cls.relative_periapsis_argument(orbit_from, orbit_to)
        theta_rel = orbit_to.theta_at_t_clock(t_clock) - orbit_from.theta_at_t_clock(t_clock) + omega_rel
        # Normalize between -180 and +180
        theta_rel = ((theta_rel + 180) % 360) - 180
        return theta_rel

    @classmethod
    def relative_theta_loop_period(cls, orbit_from, orbit_to, from_t_clock=0, max_iter=1000, tol=1e-6):
        """
        Finds the period of the relative theta loop between two coplanar (or nearly coplanar) orbits.
        """
        theta_rel0 = cls.relative_theta_at_t_clock(orbit_from, orbit_to, from_t_clock)

        def objective(t_clock):
            theta_rel = cls.relative_theta_at_t_clock(orbit_from, orbit_to, t_clock)
            diff = ((theta_rel - theta_rel0 + 180) % 360) - 180
            return abs(diff)
        
        guess = from_t_clock + abs(orbit_from.T * orbit_to.T / (orbit_to.T - orbit_from.T))
        result = minimize_scalar(
            objective,
            bracket=(guess - orbit_from.T/2, guess, guess + orbit_from.T/2),
            method='brent',
            options={'xtol': tol, 'maxiter': max_iter}
        )
        if result.success:
            return result.x
        else:
            raise RuntimeError(f"Failed to find a solution that satisfies the desired relative angle. Best difference found: {result.fun} degrees")

    @classmethod
    def find_t_clock_for_relative_theta(cls, orbit_from, orbit_to, theta_rel, from_t_clock=0, max_iter=1000, tol=1e-6):
        """
        Finds the least time t_clock for which the relative theta between two orbits is equal to theta_rel.
        Theta_rel is the relative angle of the orbit_to with respect to the orbit_from, positive if ahead, negative if behind.
        """
        def objective(t_clock):
            theta_rel_calc = cls.relative_theta_at_t_clock(orbit_from, orbit_to, t_clock)
            diff = ((theta_rel_calc - theta_rel + 180) % 360) - 180
            return abs(diff)

        loop_period = cls.relative_theta_loop_period(orbit_from, orbit_to, from_t_clock)
        
        result = minimize_scalar(
            objective,
            bounds=(from_t_clock, from_t_clock + loop_period),
            method='bounded',
            options={'xatol': tol, 'maxiter': max_iter}
        )
        
        if result.success:
            return result.x
        else:
            raise RuntimeError(f"Failed to find a solution that satisfies the desired relative angle. Best difference found: {result.fun} degrees")
        
    @staticmethod
    def delta_v_from_v_vecs(v_vec_from, v_vec_to):
        return np.linalg.norm(v_vec_to - v_vec_from)
