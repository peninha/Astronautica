# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:19:43 2020

@author: Pena
"""

from orbital_mechanics import *


""" Orbital Parameters """

Probe = {
    "father": Earth,
    "orbit": {}
    }

Moon["orbit"]["l"] = Moon["orbit"]["a"]*(1 - Moon["orbit"]["e"]**2)
Moon["orbit"]["h"] = math.sqrt(Earth["mu"] * Moon["orbit"]["l"])

RMoon, ThetaMoon, AltMoon, VMoon = [
    {} for n in range(4)]
RProbe, ThetaProbe, AltProbe, VProbe, EProbe, Mprobe = [
    {} for n in range(6)]


""" Inputs """

AltMoon["t0"] = 394091000  # Moon Altitute in meters
dRdTheta_sign = 1  # Altitude is increasing with time (+1) or Decreasing (-1)?
AltProbe["t0"] = 199970 # 200km
VProbe["t0"] = 11100 # m/s

print("Check if the Moon Velocity is: {:,.2f}".format(VFromAlt(AltMoon["t0"],
                                                               Moon)))


""" Setup """

OrbitFromPeriVandAlt(VProbe["t0"], AltProbe["t0"], Probe)
RProbe["t0"] = RFromAlt(AltProbe["t0"], Probe)
RProbe["t1"] = Moon["orbit"]["a"]
ThetaMoon0 = ThetaFromAlt(AltMoon["t0"], Moon, dRdTheta_sign)
GammaMoon0 = GammaFromTheta(ThetaMoon0, Moon)

def Equations(vector):  # Problem to Solve
    t1 = vector[0]
    Gamma1 = vector[1]
    Probe["orbit"]["argPeriapsis"] = vector[2]
    R = vector[3]

    f1 = (ThetaFromTime(t1+t0, Moon) + Moon["orbit"]["argPeriapsis"]
          )%360 - Gamma1
    f2 = (ThetaFromTime(t1, Probe) + Probe["orbit"]["argPeriapsis"]
          )%360 - Gamma1
    f3 = RFromTime(t1+t0, Moon) - R
    f4 = RFromTime(t1, Probe) - R
    return np.array([f1, f2, f3, f4])


""" Initial Guesses """

t1 = TimeFromR(Moon["orbit"]["a"], Probe)
t0 = TimeFromTheta(ThetaFromAlt(AltMoon["t0"], Moon, dRdTheta_sign), Moon)
Gamma1 = ( GammaFromTheta(ThetaFromTime(t0+t1, Moon), Moon) )%360
argPeri = ( Gamma1 - 180 )%360
R = Moon["orbit"]["a"]

print("Initial Guesses:")
print("t1: {}; Gamma1: {}; argPeri: {}; R: {}".format(
    t1, Gamma1, argPeri, R))

vector = np.array([t1, Gamma1, argPeri, R])


""" Solving """

solution = root(Equations, vector)
phaseAngle = (GammaMoon0 - Probe["orbit"]["argPeriapsis"])%360


""" Printing Results """

print("Flight Time: {:,.0f} s | {:,.1f} h | {:,.1f} d".format(
    solution.x[0], solution.x[0]/3600, solution.x[0]/(24*3600)))
print("Impact Gamma Angle: {:0.2f} degrees".format(solution.x[1]))
print("Argument of Prob's Periapsis: {:0.2f} degrees".format(solution.x[2]))
print("Altitude of Impact: {:,.0f} m".format(AltFromR(solution.x[3], Moon)))
print("Phase Angle for Burn {:0.2f} degrees".format(phaseAngle))
print("Impact Velocity: {:0.2f} m/s".format(VFromR(solution.x[3], Probe)))













"""
import math
import numpy as np
from scipy.optimize import root




def ThetaFromR(r, body, drdtheta_sign = 1):
    return drdtheta_sign * math.degrees(math.acos(
        (body["orbit"]["l"]/r - 1) / body["orbit"]["e"] )) % 360

def ThetaFromAlt(alt, body, drdtheta_sign = 1):
    return ThetaFromR(RFromAlt(alt), body, drdtheta_sign)

def ThetaFromE(E, body):
    if body["orbit"]["e"] > 1:
        return math.degrees( 2*math.atan(
        math.sqrt( abs((1 + body["orbit"]["e"]) / (1 - body["orbit"]["e"])) )
              *math.tanh(E/2)
        ) ) % 360
    else:
        return math.degrees( 2*math.atan(
        math.sqrt( abs((1 + body["orbit"]["e"]) / (1 - body["orbit"]["e"])) )
              *math.tan(E/2)
        ) ) % 360

def ThetaFromTime(t, body):
    M = MFromTime(t, body)
    E = EFromM(M, body)
    return ThetaFromE(E, body)

def ThetaFromGamma(gamma, body):
    return (gamma - body["orbit"]["argPeriapsis"])%360

def RFromTheta(theta, body):
    return body["orbit"]["l"]/(1 + body["orbit"]["e"]
                               *math.cos(math.radians(theta)))

def RFromTime(t, body):
    return RFromTheta( ThetaFromTime(t, body), body )

def RFromAlt(alt):
    return alt + Earth["radius"]

def VFromR(r, body):
    return math.sqrt( Earth["mu"] * (2/r - 1/body["orbit"]["a"]))

def VFromAlt(alt, body):
    return VFromR(RFromAlt(alt), body)

def AltFromR(r):
    return r - Earth["radius"]

def MFromTime(t, body):
    return math.sqrt( abs(Earth["mu"]/body["orbit"]["a"]**3) )*t

def MFromE(E, body):
    if body["orbit"]["e"] > 1:
        return body["orbit"]["e"]*math.sinh(E) - E
    else:
        return E - body["orbit"]["e"]*math.sin(E)

def EFromM(M, body):
    x0 = M
    result = root(KeplerEquation, x0, args=(M, body))
    return result.x[0]

def EFromTheta(theta, body):
    x = math.sqrt( abs((1 - body["orbit"]["e"])
                       / (1 + body["orbit"]["e"]))
                  ) * math.tan(math.radians(theta)/2)
    if body["orbit"]["e"] > 1:
        return math.atanh(x) * 2
    else:
        return math.atan(x) * 2

def TimeFromM(M, body):
    time = M * math.sqrt( abs(body["orbit"]["a"]**3 / Earth["mu"]) )
    if time < 0:
        return time + OrbitPeriod(body)
    else:
        return time

def TimeFromTheta(theta, body):
    E = EFromTheta(theta, body)
    M = MFromE(E, body)
    return TimeFromM(M, body)

def TimeFromR(R, body):
    return TimeFromTheta( ThetaFromR(R, body), body )

def GammaFromTheta(theta, body):
    return (theta + body["orbit"]["argPeriapsis"])%360

def KeplerEquation(E, M, body):
    if body["orbit"]["e"] > 1:
        return body["orbit"]["e"]*math.sinh(E) - E - M
    else:
        return E - body["orbit"]["e"]*math.sin(E) - M

def OrbitPeriod(body):
    return 2*math.pi * math.sqrt( abs(body["orbit"]["a"]**3 / Earth["mu"]) )

def OrbitFromPeriVandAlt(periV, periAlt, body):
    rp = RFromAlt(periAlt)
    a = 1 / (2/rp - periV**2/Earth["mu"])
    e = 1 - rp/a
    l = a * (1 - e**2)
    h = math.sqrt( Earth["mu"] * l )
    body["orbit"]["a"] = a
    body["orbit"]["e"] = e
    body["orbit"]["l"] = l
    body["orbit"]["h"] = h
    return {"a": a, "e": e, "l": l, "h": h}



G = 6.6743E-11

Sun = {
    "mu": 1.32712e+20,
    "radius": 696342000,
    }

Earth = {
    "father": Sun,
    "mu": 3.986004418e+14,
    "radius": 6371000,
    "rotationPeriod": 86164.098903691,
    "orbit": {
        "a": 149598261150.4425,
        "e": 0.01609636160505683,
        "inclination": 23.44603795469773,
        "meanAnomalyEpoch": 357.0607464120944,
        "longAscNode": 359.9965004168758,
        "argPeriapsis": 102.9720683296131,
        }
    }

Moon = {
    "father": Earth,
    "radius": 1737100,
    "mass": 7.34767309e+22,
    "rotationPeriod": 2360584.68479999,
    "orbit": {
        "a": 384308437.7707066,
        "e": 0.05328149353682574,
        "inclination": 28.36267790798491,
        "meanAnomalyEpoch": 222.7012350930954,
        "longAscNode": 2.296616161126016,
        "argPeriapsis": 199.7640930160823,
        }
    }
"""