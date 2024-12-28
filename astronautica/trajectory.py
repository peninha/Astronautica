from .orbit import Orbit

class Trajectory:
    """
    Class representing a trajectory.
    """
    def __init__(self, orbit0, t0_clock=None, theta0=None, name="Initial Orbit"):
        if not isinstance(orbit0, Orbit):
            raise ValueError("Orbit0 must be an instance of the Orbit class")
        
        self.orbits = []
        self.add_orbit(orbit0, 0, name)
        self.add_trajectory_position(0, t0_clock, theta0, name="Initial Position")

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

    def add_orbit(self, orbit, orbit_number, name="Orbit"):
        """
        Adds a new orbit to the trajectory.
        """
        if not isinstance(orbit, Orbit):
            raise ValueError("orbit must be an instance of the Orbit class")
            
        # Convert self.orbits to a dictionary if it is not already
        if isinstance(self.orbits, list):
            self.orbits = {}
            
        self.orbits[orbit_number] = {
            "orbit": orbit,
            "name": name,
            "trajectory_positions": []  # List to store trajectory positions
        }

    def add_trajectory_position(self, orbit_number=0, t_clock=None, theta=None, name=None):
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
            "name": name
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