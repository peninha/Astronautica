parameter a to ship.
parameter b to target.

local file to "orbit.txt".
deletepath(file).

function print_and_log{
    parameter string.

    print string.
    log string to file.
}

CLEARSCREEN.

print_and_log("--------------------------------------").
print_and_log("-  " + a:name).
print_and_log("--------------------------------------").

print_and_log("       Periapsis: " + round(a:orbit:periapsis/1000, 5) + " km").
print_and_log("        Apoapsis: " + round(a:orbit:apoapsis/1000, 5) + " km").
print_and_log("    Eccentricity: " + round(a:orbit:eccentricity, 5)).
print_and_log(" Semi-major axis: " + round(a:orbit:semimajoraxis/1000, 5) + " km").
print_and_log("     Inclination: " + round(a:orbit:inclination, 5) + " °").
print_and_log("             LAN: " + round(a:orbit:lan, 5) + " °").
print_and_log("Argument of peri: " + round(a:orbit:argumentofperiapsis, 5) + " °").
print_and_log("    True anomaly: " + round(a:orbit:trueanomaly, 5) + " °").
print_and_log("          Period: " + round(a:orbit:period, 0) + " s").

print_and_log(" ").

print_and_log("--------------------------------------").
print_and_log("-  " + b:name).
print_and_log("--------------------------------------").

print_and_log("       Periapsis: " + round(b:orbit:periapsis/1000, 5) + " km").
print_and_log("        Apoapsis: " + round(b:orbit:apoapsis/1000, 5) + " km").
print_and_log("    Eccentricity: " + round(b:orbit:eccentricity, 5)).
print_and_log(" Semi-major axis: " + round(b:orbit:semimajoraxis/1000, 5) + " km").
print_and_log("     Inclination: " + round(b:orbit:inclination, 5) + " °").
print_and_log("             LAN: " + round(b:orbit:lan, 5) + " °").
print_and_log("Argument of peri: " + round(b:orbit:argumentofperiapsis, 5) + " °").
print_and_log("    True anomaly: " + round(b:orbit:trueanomaly, 5) + " °").
print_and_log("          Period: " + round(b:orbit:period, 0) + " s").

print_and_log(" ").

// but does it by total seconds (2*60 + 30 = 150):
//local maneuver to NODE(time:seconds + 150, 0, 50, 10).
//add maneuver.

kuniverse:pause().