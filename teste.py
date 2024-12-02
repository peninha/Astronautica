from orbitalmechanics import Orbit

M_earth = 5.9722e24 # [kg]
rp = 7000 # [km]    
e = 0 # [km]

try:
    orbit_rp_e = Orbit(m1=M_earth, rp=rp, e=e)
    print(orbit_rp_e)

    orbit_rp_ra = Orbit(m1=M_earth, rp=rp, ra=orbit_rp_e.ra)
    print(orbit_rp_ra)

    orbit_rp_h = Orbit(m1=M_earth, rp=rp, h=orbit_rp_e.h)
    print(orbit_rp_h)

    orbit_e_h = Orbit(m1=M_earth, e=orbit_rp_e.e, h=orbit_rp_e.h)
    print(orbit_e_h)

    orbit_a_e = Orbit(m1=M_earth, a=orbit_rp_e.a, e=orbit_rp_e.e)
    print(orbit_a_e)

except ValueError as e:
    print(e)
