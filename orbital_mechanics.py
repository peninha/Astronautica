# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:19:43 2020

@author: Pena
"""
import math
import numpy as np
from scipy.optimize import root


""" Orbit Functions """

def ThetaFromR(r, body, drdtheta_sign = 1):
    return drdtheta_sign * math.degrees(math.acos(
        (body["orbit"]["l"]/r - 1) / body["orbit"]["e"] )) % 360

def ThetaFromAlt(alt, body, drdtheta_sign = 1):
    return ThetaFromR(RFromAlt(alt, body), body, drdtheta_sign)

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

def RFromAlt(alt, body):
    return alt + body["father"]["radius"]

def VFromR(r, body):
    return math.sqrt( body["father"]["mu"] * (2/r - 1/body["orbit"]["a"]))

def VFromAlt(alt, body):
    return VFromR(RFromAlt(alt, body), body)

def AltFromR(r, body):
    return r - body["father"]["radius"]

def MFromTime(t, body):
    #return math.sqrt( abs(body["father"]["mu"]/body["orbit"]["a"]**3) )*t
    return t / OrbitPeriod(body) * 2*math.pi

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
    time = M / (2*math.pi) * OrbitPeriod(body)
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
    if "period" in body["orbit"]:
        return body["orbit"]["period"]
    else:
        return 2*math.pi * math.sqrt( abs(body["orbit"]["a"]**3
                                      / body["father"]["mu"]) )

def OrbitFromPeriVandAlt(periV, periAlt, body):
    rp = RFromAlt(periAlt, body)
    a = 1 / (2/rp - periV**2/body["father"]["mu"])
    e = 1 - rp/a
    l = a * (1 - e**2)
    h = math.sqrt( body["father"]["mu"] * l )
    body["orbit"]["a"] = a
    body["orbit"]["e"] = e
    body["orbit"]["l"] = l
    body["orbit"]["h"] = h
    return {"a": a, "e": e, "l": l, "h": h}


""" Orbital Parameters """
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
        #"period": 2356544.495,
        #"period": 2360584.68479999,
        }
    }

Mars = {
    "father": Sun,
    "radius": 3375800,
    "mu": 4.282831e+13,
    "rotationPeriod": 88642.6848,
    "orbit": {
        "a": 227949699961.9763,
        "e": 0.09326110278323557,
        "inclination": 24.69272426910055,
        "meanAnomalyEpoch": 169.3913127942378,
        "longAscNode": 3.351911063089117,
        "argPeriapsis": 332.1022655295414,
        }
    }

Bodies = [Earth, Mars, Moon]

for body in Bodies:
    if "mu" not in body:
        body["mu"] = body["mass"] * G
    if "mass" not in body:
        body["mass"] = body["mu"] / G
    if "SOI" not in body:
        body["SOI"] = body["orbit"]["a"]*(body["mu"]
                                          / body["father"]["mu"])**(2/5)
    #body["orbit"]["period"] = 2*math.pi * math.sqrt(body["orbit"]["a"]**3
     #                                     / body["father"]["mu"])
    body["orbit"]["apoapsis"] = AltFromR(body["orbit"]["a"]
                                         *(1 + body["orbit"]["e"]), body)
    body["orbit"]["periapsis"] = AltFromR(body["orbit"]["a"]
                                         *(1 - body["orbit"]["e"]), body)
    body["orbit"]["l"] = body["orbit"]["a"]*(1 - body["orbit"]["e"]**2)
    body["orbit"]["h"] = math.sqrt(body["father"]["mu"] * body["orbit"]["l"])
