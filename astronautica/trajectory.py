from .orbit import Orbit
from .maneuver import Maneuver

class Trajectory:
    """
    Class representing a trajectory.
    """
    def __init__(self, orbit0, orbit_name="Initial Orbit", position_name="Initial position"):
        if not isinstance(orbit0, Orbit):
            raise ValueError("Orbit0 must be an instance of the Orbit class")
        
        self.orbits = []
        self.add_orbit(orbit0, 0, orbit_name)
        self.add_trajectory_position(0, orbit0.t0_clock, orbit0.theta0, name=position_name)

    def __str__(self):
        # Get all attributes of the class
        attributes = vars(self)
        # Initial string
        result = "Trajectory Parameters:\n"
        # Add each attribute to the string
        for name, value in attributes.items():
            # Format numbers as float with 4 decimal places
            if isinstance(value, (int, float)):
                result += f"{name}: {value:.4f}\n"
            else:
                result += f"{name}: {value}\n"

        return result

    def add_orbit(self, orbit, orbit_number, name=None):
        """
        Adds a new orbit to the trajectory.
        """
        if not isinstance(orbit, Orbit):
            raise ValueError("orbit must be an instance of the Orbit class")
        
        if name is not None:
            orbit.name = name
            
        # Convert self.orbits to a dictionary if it is not already
        if isinstance(self.orbits, list):
            self.orbits = {}
            
        self.orbits[orbit_number] = {
            "orbit": orbit,
            "trajectory_positions": [],  # List to store trajectory positions
            "maneuvers": []  # List to store maneuvers
        }

    def add_maneuver(self, orbit_number, maneuver, name=None, orbit_name=None, position_name=None):
        if not isinstance(maneuver, Maneuver):
            raise ValueError("maneuver must be an instance of the Maneuver class")
        
        if name is None:
            name = maneuver.name
        if orbit_name is None:
            orbit_name = maneuver.orbit_name
        if position_name is None:
            position_name = maneuver.position_name
        
        self.orbits[orbit_number]["maneuvers"].append(maneuver)
        orbit = self.orbits[orbit_number]["orbit"]
        theta_node = orbit.theta_at_t_clock(maneuver.t_clock)
        self.add_trajectory_position(orbit_number, maneuver.t_clock, theta_node, name=position_name)

        if maneuver.mode == "RTN":
            delta_v_vec_bc = orbit.convert_RTN_to_cartesian(maneuver.RTN, maneuver.t_clock, frame="bodycentric")
        elif maneuver.mode == "RPN":
            delta_v_vec_bc = orbit.convert_RPN_to_cartesian(maneuver.RPN, maneuver.t_clock, frame="bodycentric")
        r_vec_bc, v_vec_bc = orbit.state_vectors_at_t_clock(maneuver.t_clock, frame="bodycentric")
        new_v_vec_bc = v_vec_bc + delta_v_vec_bc
        
        maneuver.delta_v_vec_bc = delta_v_vec_bc
        new_orbit = Orbit.from_state_vectors(orbit.main_body, r_vec_bc, new_v_vec_bc, t0_clock=maneuver.t_clock, name=orbit_name, position_name=position_name)
        self.add_orbit(new_orbit, orbit_number + 1, name=orbit_name)
        self.add_trajectory_position(orbit_number + 1, maneuver.t_clock, theta_node, name=position_name)
        return new_orbit

    def add_trajectory_position(self, orbit_number=0, t_clock=None, theta=None, name="Position"):
        """
        Adds a trajectory position to a specific orbit in the trajectory.
        """
        if orbit_number not in self.orbits:
            raise ValueError(f"Orbit number {orbit_number} not found")
            
        if t_clock is None and theta is None:
            raise ValueError("t_clock and theta cannot be None simultaneously")
            
        orbit = self.orbits[orbit_number]["orbit"]
        
        if t_clock is None:
            t_clock = orbit.t_clock_at_theta(theta)
        if theta is None:
            theta = orbit.theta_at_t_clock(t_clock)

        position = {
            "t_clock": t_clock,
            "theta": theta,
            "name": str(orbit_number) + " - " + name
        }
        
        self.orbits[orbit_number]["trajectory_positions"].append(position)
        self.sort_trajectory_positions(orbit_number)

    def sort_trajectory_positions(self, orbit_number):
        self.orbits[orbit_number]["trajectory_positions"].sort(key=lambda x: x["t_clock"])

    def get_trajectory_position(self, orbit_number, position_index="last"):
        if position_index == "last":
            position_index = -1
        elif position_index == "first":
            position_index = 0
        elif isinstance(position_index, int):
            if position_index < 0 or position_index >= len(self.orbits[orbit_number]["trajectory_positions"]):
                raise ValueError("Invalid position index")
        else:
            raise ValueError("Invalid position index")
        return self.orbits[orbit_number]["trajectory_positions"][position_index]
    
    def state_vectors_at_trajectory_position(self, orbit_number, position_index="last", frame="bodycentric"):
        position = self.get_trajectory_position(orbit_number, position_index)
        return self.orbits[orbit_number]["orbit"].state_vectors_at_t_clock(position["t_clock"], frame)