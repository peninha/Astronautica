from .body import Body
import numpy as np

class Groundstations:
    """
    Class representing groundtrack stations on the surface of a body.
    """
    def __init__(self, body = Body("earth")):
        """
        Initialize the main body.
        """
        self.body = body
        self.stations = []

    def add_groundstation(self, latitude = 0.0, longitude = 0.0, name = "Groundstation"):
        """
        Add a groundstation to the list.
        """
        self.stations.append({"latitude": latitude, "longitude": longitude, "name": name})

