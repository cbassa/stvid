#!/usr/bin/env python3
import numpy as np
from scipy.interpolate import interpolate, CubicSpline

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, get_body, AltAz, ITRS, GCRS


def compute_reference_distance(t, loc, az, alt, height):
    # Distance to search over
    distance = np.linspace(0, 5 * height.value) * height.unit
    asat = AltAz(obstime=t, location=loc, az=az, alt=alt, distance=distance)
    itrs = asat.transform_to(ITRS(obstime=t))
    l = EarthLocation.from_geocentric(*itrs.cartesian.xyz)
    fint = interpolate.interp1d(l.height.to(height.unit), distance.to(height.unit))

    return fint(height) * height.unit


def spline_roots(x, y, ydescend, yascend):
    for y0, sign in zip([ydescend, yascend], [-1, +1]):
        # Define spline
        func = CubicSpline(x, y - y0)

        # Find roots
        xr = func.roots(extrapolate=False)

        # Select based on sign
        c = np.sign(func.derivative()(xr)) == sign

        if sign == +1:
            xascend = xr[c]
        elif sign == -1:
            xdescend = xr[c]

    return xascend, xdescend


def shadow_angle(t, loc, gsun, az, alt, height):
    # Distance
    dist = compute_reference_distance(t[0], loc, az, alt, height)

    # Compute satellite position
    gsat = AltAz(obstime=t, location=loc, az=az, alt=alt, distance=dist).transform_to(
        GCRS(obstime=t)
    )
    psat = gsat.cartesian
    psun = gsun.cartesian
    dp = psun - psat

    # Compute satellite offset from horizon
    angle = np.arccos(-dp.dot(psat) / (dp.norm() * psat.norm()))
    angle_earth = np.arcsin(6378.135 * u.km / psat.norm())

    return angle - angle_earth


def observe_logic(tnow, loc, sunaltset, sunaltrise, az, alt, height):
    # Next 48h to ensure roots are found
    t = tnow + np.linspace(0, 2, 64) * u.d

    # Sun location
    gsun = get_body("sun", t, location=loc)
    asun = gsun.transform_to(AltAz(obstime=t, location=loc))

    # Parameters

    # Solar altitude constraints
    sunalt = asun.alt
    if np.min(sunalt) > min(sunaltset, sunaltrise):
        sunstate = "Sun never sets"
    elif np.max(sunalt) < max(sunaltset, sunaltrise):
        sunstate = "Sun never rises"
    else:
        sunstate = "Sun rises and sets"
    trise, tset = Time(spline_roots(t.mjd, sunalt, sunaltset, sunaltrise), format="mjd")

    # Aimpoint constraints
    if az != None and alt != None and height != None:
        angle = shadow_angle(t, loc, gsun, az, alt, height)
        if np.min(angle) > 0:
            shadowstate = "Aimpoint is always sunlit"
        elif np.max(angle) < 0:
            shadowstate = "Aimpoint is never sunlit"
        else:
            shadowstate = "Aimpoint eclipses"
        tegress, tingress = spline_roots(t.mjd, angle, 0, 0)
        if len(tegress) > 0:
            tegress = Time(tegress, format="mjd")
        if len(tingress) > 0:
            tingress = Time(tingress, format="mjd")
    else:
        shadowstate = "Aimpoint is always sunlit"
        tegress, tingress = [], []

    # Get first value
    if sunstate == "Sun rises and sets":
        trise, tset = np.sort(trise)[0], np.sort(tset)[0]
    if shadowstate == "Aimpoint eclipses":
        tingress, tegress = np.sort(tingress)[0], np.sort(tegress)[0]

    # Logic
    if sunstate == "Sun rises and sets" and shadowstate == "Aimpoint eclipses":
        if tset < tingress < tegress < trise:
            state = f"The Sun is above the horizon. Sunset at {tset.isot}."
            wait_time = np.floor((tset - tnow).to(u.s).value)
            action = "wait"
            tend = tingress
        elif tingress < tegress < trise < tset:
            state = "The Sun is below the horizon and the aimpoint is sunlit."
            wait_time = 0
            action = "observe"
            tend = tingress
        elif tegress < trise < tset < tingress:
            state = f"The aimpoint is not sunlit. Shadow egress at {tegress.isot}"
            wait_time = np.floor((tegress - tnow).to(u.s).value)
            action = "wait"
            tend = trise
        elif trise < tset < tingress < tegress:
            state = "The Sun is below the horizon and the aimpoint is sunlit."
            wait_time = 0
            action = "observe"
            tend = trise
    elif (sunstate == "Sun rises and sets" and shadowstate == "Aimpoint is always sunlit"):
        if tset < trise:
            state = f"The Sun is above the horizon. Sunset at {tset.isot}."
            wait_time = np.floor((tset - tnow).to(u.s).value)
            action = "wait"
            tend = trise
        elif trise < tset:
            state = "The Sun is below the horizon."
            wait_time = 0
            action = "observe"
            tend = trise
    elif sunstate == "Sun never rises" and shadowstate == "Aimpoint eclipses":
        if tingress < tegress:
            state = "The Sun is below the horizon and the aimpoint is sunlit."
            wait_time = 0
            action = "observe"
            tend = tingress
        elif tegress < tingress:
            state = "The Sun is below the horizon and the aimpoint is not sunlit."
            wait_time = np.floor((tegress - tnow).to(u.s).value)
            action = "wait"
            tend = trise
    elif sunstate == "Sun never sets":
        state = "The Sun never sets."
        wait_time = 86400
        action = "wait"
        tend = tnow + 24 * u.h
    elif shadowstate == "Aimpoint is never sunlit":
        state = "The aimpoint is never sunlit."
        wait_time = 86400
        action = "wait"
        tend = tnow + 24 * u.h
    else:
        state = "Unknown state."
        wait_time = 86400
        action = "wait"
        tend = tnow + 24 * u.h

    return action, wait_time, tend, state
