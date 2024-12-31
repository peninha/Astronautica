from astronautica import Orbit, Body, Plotter, Maneuver
import numpy as np

"""
Spacecraft at A and B are in the same orbit (1). At the instant shown in Fig. 6.11 the chaser vehicle
at A executes a phasing maneuver so as to catch the target spacecraft back at A after just one revolution
of the chaser’s phasing orbit (2). What is the required total delta-v?
"""
earth = Body("earth")

rp = 6800
ra = 13600

orbita = Orbit.from_apsis(earth, rp, ra)

rp = orbita.rp

delta_v1, delta_v2 = Maneuver.delta_v_for_phase_change(orbita, 90, theta_burn=0, n=1)

print(f"Delta-v 1: {delta_v1} km/s")
print(f"Delta-v 2: {delta_v2} km/s")
print(f"Total delta-v: {np.abs(delta_v1) + np.abs(delta_v2)} km/s")