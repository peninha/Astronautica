from astronautica import Orbit, Body, Plotter, Maneuver, Trajectory
import numpy as np

"""
A geocentric satellite in orbit 1 of Fig. 6.15 executes a delta-v maneuver at A, which places it on orbit 2, for reentry at D.
Calculate Δv at A and its direction relative to the local horizon
"""
earth = Body("earth")
