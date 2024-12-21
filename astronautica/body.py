# astronautica/body.py
class Body:
    """
    Class representing a celestial body.
    It can be initialized with a predefined name (Earth, Moon, etc.) or with custom parameters.
    """
    _known_bodies = {
        "spherical_earth": {"mass": 5.972e24, "radius": 6378.137, "rotation_speed": 7.292115e-5, "J2": 0, "flatness": 0},
        "earth": {"mass": 5.972e24, "radius": 6378.137, "rotation_speed": 7.292115e-5, "J2": 1.0826e-3, "flatness": 0.0033528},
        "moon": {"mass": 7.34767309e22, "radius": 1737.4, "rotation_speed": 2.6617e-6, "J2": 0.0, "flatness": 0.0},
        "mars": {"mass": 6.4171e23, "radius": 3389.5, "rotation_speed": 7.088e-5, "J2": 0.00589, "flatness": 0.00589},
    }

    def __init__(self, name=None, mass=0, radius=0, rotation_speed=0, J2=0, flatness=0):
        if name in self._known_bodies:
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

    def __repr__(self):
        return (
            f"Body(name='{self.name}', mass={self.mass} kg, radius={self.radius} km, "
            f"rotation_speed={self.rotation_speed} rad/s, J2={self.J2}, flatness={self.flatness})"
        )

# Example usage
if __name__ == "__main__":
    earth = Body(name="Earth")
    print(earth)

    moon = Body(name="Moon")
    print(moon)

    custom_body = Body(name="Exoplanet X", mass=6e24, radius=7000, rotation_speed=1e-5, J2=0.002, flatness=0.003)
    print(custom_body)
