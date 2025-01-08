import numpy as np

class Body:
    """
    Class representing a celestial body.
    It can be initialized with a predefined name (Earth, Moon, etc.) or with custom parameters.
    """

    ########### CONSTANTS ###########
    G = 6.67430e-20  # Gravitational constant in km^3/kg/s^2
 
    ########### KNOWN BODIES ###########
    _known_bodies = {
        "spherical_earth": {"mass": 5.972e24, "radius": 6378.137, "rotation_speed": 0.004178074216, "J2": 0, "flatness": 0},
        "earth": {"mass": 5.972e24,
                   "radius": 6378.137,
                   #"radius": 6378,
                   "rotation_speed": 0.004178074216,
                   "J2": 1.08262668e-3,
                   #"J2": 1.08263e-3,
                   "flatness": 0.0033528},
        "squished_earth": {"mass": 5.972e24,
                   "radius": 6378.137,
                   #"radius": 6378,
                   "rotation_speed": 0.004178074216,
                   "J2": 1.08262668e-3 * 100,
                   #"J2": 1.08263e-3,
                   "flatness": 0.0033528 * 100},
        "moon": {"mass": 7.34767309e22, "radius": 1737.4, "rotation_speed": 0.0001525041, "J2": 0.0, "flatness": 0.0},
        "mars": {"mass": 6.4171e23, "radius": 3389.5, "rotation_speed": 0.004061248134, "J2": 0.00589, "flatness": 0.00589},
    }

    def __init__(self, name=None, mass=0, radius=0, rotation_speed=0, J2=0, flatness=0):
        if isinstance(name, str) and name in self._known_bodies:
            self.name = name
            self.mass = self._known_bodies[name]["mass"]
            self.radius = self._known_bodies[name]["radius"]
            self.rotation_speed = self._known_bodies[name]["rotation_speed"]
            self.J2 = self._known_bodies[name]["J2"]
            self.flatness = self._known_bodies[name]["flatness"]
        else:
            self.name = name if name else "Custom Body"
            self.mass = mass
            self.radius = radius
            self.rotation_speed = rotation_speed
            self.J2 = J2
            self.flatness = flatness
   
    def __str__(self):
        # Get all attributes of the class
        attributes = vars(self)
        # Initial string
        result = "Body Parameters:\n"
        # Add each attribute to the string
        for name, value in attributes.items():
            # Format numbers as float with 4 decimal places
            if isinstance(value, (int, float)):
                result += f"{name}: {value:.4f}\n"
            else:
                result += f"{name}: {value}\n"

        return result

    def altitude(self, r):
        """
        Calculate the altitude of a point above the body's surface.
        """
        return r - self.radius
    
    def radius_from_altitude(self, altitude):
        """
        Calculate the radius of a point above the body's surface.
        """
        return self.radius + altitude

    def v_escape(self, r):
        """
        Calculate the escape velocity from a point above the body's surface.
        """
        return np.sqrt(2 * self.G * self.mass / r)

    def mu(self):
        """
        Calculate the gravitational parameter of the body.
        """
        return 398600
        #return self.G * self.mass

    def geosynchronous_radius(self):
        """
        Calculate the geosynchronous radius of the body.
        """
        T = 360 / self.rotation_speed
        return (self.mu() * T**2 / (4 * np.pi**2))**(1/3)
    
    def geosynchronous_altitude(self):
        """
        Calculate the geosynchronous altitude of the body.
        """
        return self.altitude(self.geosynchronous_radius())

