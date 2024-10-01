# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:19:43 2020

@author: Pena
"""

from orbital_mechanics import *

def YearDays(m, d):
    months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    days = 0
    for month in range(m-1):
        days += months[month]
    days += d
    return days

def TimeFromDate(Y, m, d, h, min, s):
    return s + 60*(min + 60*(h + 24*(YearDays(m, d)+2 + 365*(Y-1951))))

def DurationFromTime(s):
    y = s / (60*60*24*365)
    d = y % 1 * 365
    h = d % 1 * 24
    m = h % 1 * 60
    s = m % 1 * 60
    #return str((y))+"y "+str((d))+"d "+str(int(h))+"h "\
    return str(int(y))+"y "+str(int(d))+"d "+str(int(h))+"h "\
        +str(int(m))+"m "+str(s)+"s"

#Y = 1960
#m = 6
#d = 28
#h = 18
#min = 20
#s = 14

#t = TimeFromDate(Y, m, d, h, min, s)
#t = 299528514.931 + 1568384
#t = -31542641.784


#t += 68384.9 + 31542641.784


body = Moon
body["orbit"]["period"] = 2360584.68479999
body["orbit"]["period"] = 2360584.68479999
tadjust = 0

print("Period: {}".format(DurationFromTime(body["orbit"]["period"])))
print("")

t = 1000

t += 31542641.784
t += tadjust
theta = ThetaFromTime(t, body) + body["orbit"]["meanAnomalyEpoch"]
radius = RFromTheta(theta, body)
altitude = AltFromR(radius, body)
velocity = VFromAlt(altitude, body)

#print("Theta:", theta)
#print("Radius:", radius)
print("Altitude: {:,.0f}".format(altitude))
print("Altref:   357.468.727")
print("Velocity:", velocity)
#print("Apoapsis: {:,}".format(body["orbit"]["apoapsis"]))
#print("Periapsis: {:,}".format(body["orbit"]["periapsis"]))
print("")

t = 100000000

t += 31542641.784
t += tadjust
theta = ThetaFromTime(t, body) + body["orbit"]["meanAnomalyEpoch"]
radius = RFromTheta(theta, body)
altitude = AltFromR(radius, body)
velocity = VFromAlt(altitude, body)

#print("Theta:", theta)
#print("Radius:", radius)
print("Altitude: {:,.0f}".format(altitude))
print("Altref:   397.059.452")
print("Velocity:", velocity)
#print("Apoapsis: {:,}".format(body["orbit"]["apoapsis"]))
#print("Periapsis: {:,}".format(body["orbit"]["periapsis"]))
print("")



body = Moon
t = -31542641.784
t = 1120

tepoch = -31542641.784
t = t - tepoch
meanAnomalyAdjust = 0

period = body["orbit"]["period"]
n = math.sqrt( (body["father"]["mu"]+body["mu"]) / body["orbit"]["a"]**3 )
period = 2*math.pi / n

#n = 2*math.pi / period
thetaM = n * t + body["orbit"]["meanAnomalyEpoch"] + meanAnomalyAdjust
thetaM = 180 * 3
thetaE = EFromM(thetaM, body)
theta = ThetaFromE(thetaE, body)
r = body["orbit"]["a"]*(1 - body["orbit"]["e"]*math.cos(math.radians(thetaE)))
r2 = RFromTheta(theta, body)
alt = AltFromR(r, body)
alt2 = AltFromR(r2, body)
print(period)
print(DurationFromTime(period))
print("ThetaM:", thetaM)
print("ThetaE:", thetaE)
print("Theta: ", theta)
print("r ", r)
print("r2", r2)
print("alt2: {:,.0f}".format(alt2))
print("alt:  {:,.0f}".format(alt))
print("peri: {:,.0f}".format(body["orbit"]["periapsis"]))
print("apo: {:,.0f}".format(body["orbit"]["apoapsis"]))
