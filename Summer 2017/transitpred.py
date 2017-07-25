#!/bin/bash
""":"
exec python $0 ${1+"$@"}
"""
#"

# This utility generates transit predictions for a given night or for
# a range of nights, for a fixed site, and for a system with parameters
# specified on the command-line, or for a system (or all systems) from
# the HN and/or HS databases.

import os, re, base64, getopt, sys, datetime, subprocess, math, numpy

from HATpipepy.LibNova.SiderealTime import apparent_sidereal_time, \
    mean_sidereal_time
from HATpipepy.LibNova.Lunar import lunar_equ_coords
from HATpipepy.LibNova.Solar import solar_equ_coords
from HATpipepy.LibNova.AngularSeparation import angular_separation
from HATpipepy.LibNova.Refraction import refraction
from HATpipepy.LibNova.Transform import hrz_from_equ, equ_from_hrz
from HATpipepy.LibNova.JulianDay import julian_day

isimportmysql=0

homedir=os.getenv("HOME")
try:
    testfile = open(homedir+'/.inspectcandidates/passwd',"r")
except IOError:
    # It doesn't exist, so prompt for the password
    from Tkinter import *
    import Pmw, subprocess, stat

def Usage():
    sys.stderr.write("Usage:\t"+sys.argv[0]+"\n")
    sys.stderr.write("\t<-d yyyy.mmdd | --start yyyy.mmdd --stop yyyy.mmdd>\n")
    sys.stderr.write("\t[-l flwo|keck|ftn|fts|lapalma|wise|piszkes|lasilla|ctio|ohp\n")
    sys.stderr.write("\t\t|mcdonald|saao\n")
    sys.stderr.write("\t\t| --lon longitude(degE) --lat lat --tz tz]\n")
    sys.stderr.write("\t<  [--obj HTRname|HATSname|HATname]]\n")
    sys.stderr.write("\t | [--allHN]  | [--allHS]  |  [--allTEP]\n")
    sys.stderr.write("\t | [--name name --ra ra --dec dec -p period --tc Tc -q q]>\n")
    sys.stderr.write("\t[--secondary] [--twialt alt] [--objalt alt]\n")
    exit(1)

# Class defining a transit event
class TransitEvent:
    def __init__(self, obj, prio, n, night, symbol, jdstart, ltstart, altstart, jdcent, ltcent, altcent, jdstop, ltstop, altstop, moonsep, futype, magv, starobj):
        self.obj = obj
        self.prio = prio
        self.n = n
        self.night = night
        self.symbol = symbol
        self.jdstart = jdstart
        self.ltstart = ltstart
        self.altstart = altstart
        self.jdcent = jdcent
        self.ltcent = ltcent
        self.altcent = altcent
        self.jdstop = jdstop
        self.ltstop = ltstop
        self.altstop = altstop
        self.moonsep = moonsep
        self.futype = futype
        self.magv = magv
        self.eventPriority = 0.0
        self.starobj = starobj
        self.ptotal = 0.0
        self.count = 0.0
    def priority():
        global shallowmagmin    # Needed to modify global copy of shallowmagmin
        global shallowdurmin
        global shallowdurmax

        P_alt = abs(math.cos(math.radians(90 - self.altcent)))
        P_mag = 10.0**(-0.02*(self.magv - shallowmagmin))
        if self.symbol == 'OIBEO':
            P_flag = 1.0
        elif self.symbol == '-IBEO' or self.symbol == 'OIBE-':
            P_flag = 0.75
        elif self.symbol == 'OIB--' or self.symbol == '--BEO':
            P_flag = 0.5
        elif self.symbol == 'OI---' or self.symbol == '-IBE-' or self.symbol == '---EO':
            P_flag = 0.1
        else:
            P_flag = 0.0
        P_objprio = 10.0/(2.0**int(self.starobj.prio))
        P_period = self.starobj.P/(self.ptotal/self.count) 
        P_dur = 0.25 + (0.5*(self.starobj.HalfDuration * 2) - shallowdurmin)/(shallowdurmax - shallowdurmin)
        
        minalt = math.cos(math.radians(90 - 30.0))
        maxalt = math.cos(math.radians(90 - 90.0))
        minmag = (10.0**(-0.02*(shallowmagmin - shallowmagmin)))
        maxmag = (10.0**(-0.02*(7.438 - shallowmagmin)))
        minper = 0.114018/(self.ptotal/self.count)
        maxper = 32.9719/(self.ptotal/self.count)
        mindur = 0.25 + ((0.5*(shallowdurmin - shallowdurmin))/(shallowdurmax - shallowdurmin))
        maxdur = 0.25 + ((0.5*(shallowdurmax - shallowdurmin))/(shallowdurmax - shallowdurmin))

        c_alt = math.log(1.25/0.75)/(math.log(maxalt/minalt))
        c_mag = math.log(1.5/0.5)/(math.log(maxmag/minmag))
        c_per = math.log(1.5/0.5)/(math.log(maxper/minper))
        c_dur = math.log(1.1/0.9)/(math.log(maxdur/mindur))

        P_alt = 0.75 * ((P_alt/minalt)**c_alt)
        P_mag = 0.5 * ((P_mag/minmag)**c_mag)
        P_period = 0.5 * ((P_period/minper)**c_per)
        P_dur = 0.9 * ((P_dur/mindur)**c_dur)

        self.eventPriority = P_alt * P_mag * P_flag * P_objprio * P_period * P_dur #NEED TO GET FLAG SOMEHOW
    # get_period function allows the tracking of ptotal and count, which are used in priority()
    # ptotal and count are incremented when the transitevents are created
    def get_period(per, count):
        self.ptotal = per
        self.count = count
    def __repr__(self):
        return repr((self.obj, self.prio, self.n, self.night, self.symbol,
                     self.jdstart, self.ltstart, self.altstart, self.jdcent,
                     self.ltcent, self.altcent, self.jdstop, self.ltstop,
                     self.altstop, self.moonsep, self.futype, self.magv, self.eventPriority, self.starobj))
    def printvalues(self):
        print "%s %s %d %s %5s %.17g %s %6.2f %.17g %s %6.2f %.17g %s %6.2f %7.2f %s %5g %5g" % (self.obj, self.prio, self.n, self.night, self.symbol, self.jdstart, self.ltstart, self.altstart, self.jdcent, self.ltcent, self.altcent, self.jdstop, self.ltstop, self.altstop, self.moonsep, self.futype, self.magv, self.eventPriority)

# Class storing the parameters of a transit candidate
class Star:
    def __init__(self, obj, prio, RA, Dec, P, Epoch, HalfDuration, futype, magv):
        self.obj = obj
        self.prio = prio
        self.RA = RA
        self.Dec = Dec
        self.P = P
        self.Epoch = Epoch
        self.HalfDuration = HalfDuration
        self.futype = futype
        self.magv = magv
    def __repr__(self):
        return repr((self.obj, self.prio, self.RA, self.Dec, self.P,
                     self.Epoch, self.HalfDuration, self.futype, self.magv))

def HA(jd, RA, sitelong):
    """ Returns the hour angle corresponding to the given julian date and
        RA """
    ha=(apparent_sidereal_time(jd)-(RA-sitelong)/15.0)%24.0
    if ha>12.0: ha-=24.0
    return ha

# Returns the apparent LST (in degrees) for a given JD and degree east longitude
# using libnova
def LST(jd, longitude):
    lst=(apparent_sidereal_time(jd)*15.0+longitude)%360.0
    return lst

# Calculate the LMST (in degrees) for a given JD and degree east longitude
# Taken from skycalc code, see also Allen's Astrophysical Quantities
#                                   Fourth Edition, section 2.4.
def LMST(jd, longitude):
    longit = -longitude/15.
    jdnoon2000jan1 = 2451545.0
    jdint = int(jd)
    jdfrac = jd - jdint
    if(jdfrac < 0.5):
        jdmid = jdint - 0.5
        ut = jdfrac + 0.5
    else:
        jdmid = jdint + 0.5
        ut = jdfrac - 0.5
    t = (jdmid - jdnoon2000jan1)/36525
    sid_g = (24110.54841+8640184.812866*t+0.093104*t*t-6.2e-6*t*t*t)/86400
    sid_int = int(sid_g)
    sid_g = sid_g - sid_int
    sid_g = sid_g + 1.0027379093*ut - longit/24.
    sid_int = int(sid_g)
    sid_g = (sid_g - sid_int) * 24.
    if(sid_g < 0.):
        sid_g = sid_g + 24.
    sid_g = sid_g*15.
    return(sid_g)

# Calculate the altitude of a given RA, DEC at latitude lat and LST lst
# All input quantities are in degrees.
def getalt(RA, Dec, lst, lat):
    #c1 = math.cos((RA - lst)*math.pi/180.)
    #c2 = math.cos(lat*math.pi/180.)*math.cos(Dec*math.pi/180.)
    #s1 = math.sin(lat*math.pi/180.)*math.sin(Dec*math.pi/180.)
    #altitude = 90.0 - 180.0*math.acos((c1*c2 + s1))/math.pi
    altitude = 90.0 - angular_separation(RA, Dec, lst, lat)
    return(altitude)

# Calculate the Hour-Angle at which a given altitude is achieved
# All quantities are in degrees
def alt_to_HA(alt, Dec, lat):
    c1 = math.cos((90.0 - alt)*math.pi/180.0)
    c2 = math.cos(lat*math.pi/180.0)*math.cos(Dec*math.pi/180.0)
    s1 = math.sin(lat*math.pi/180.0)*math.sin(Dec*math.pi/180.0)
    if c2 == 0:
        # Either we are at a pole, or we are looking at a pole
        # The altitude doesn't change with time in this case.
        return -1000
    val = (c1 - s1)/c2
    if val < -1.0 or val > 1.0:
        # The Dec never crosses this altitude
        return -1000
    ha = math.acos(val)*180.0/math.pi
    return ha

# Returns approximate RA and Dec of the Sun for a given JD
# RA and Dec both in degrees.
# See Astronomical Almanac, also skycalc code
def SunCoords(jd):
    jdnoon2000jan1 = 2451545.0
    n = jd - jdnoon2000jan1
    L = 280.460 + 0.9856474*n
    g = 357.528 + 0.9856003*n
    if(L < 0 or L >= 360.0):
        L = L - 360.0*int(L/360.0)
    if(g < 0 or g >= 360.0):
        g = g - 360.0*int(g/360.0)
    lam = L + 1.915*math.sin(g*math.pi/180.0) + 0.20*math.sin(g*2.0*math.pi/180.0)
    eps = 23.439 - 0.0000004*n
    x = math.cos(lam*math.pi/180.0)
    y = math.cos(eps*math.pi/180.0)*math.sin(lam*math.pi/180.0)
    z = math.sin(eps*math.pi/180.0)*math.sin(lam*math.pi/180.0)
    raval = 180.0*math.atan2(y, x)/math.pi
    if(raval < 0.):
        raval = raval + 360.
    decval = math.asin(z)*180.0/math.pi
    return [raval, decval]

# Determine the JD at which the Sun is at a given altitude given an 
# initial guess (see skycalc.v5.c)
# Alt, lat and longitude in degrees. Longit is degrees west.
def jd_sun_alt(alt, jdin, lat, longit):
    delta = 0.002
    alterr = 0.001
    maxstep = 10000
    i = 0
    
    jdguess = jdin
    
    #[raval, decval] = SunCoords(jdguess)
    [raval, decval] = solar_equ_coords(jdguess)
    lstval = LST(jdguess,longit)
    alt2 = getalt(raval, decval, lstval, lat)
    if abs(alt2 - alt) < alterr:
        return jdin
    jdguess = jdguess + delta
    #[raval, decval] = SunCoords(jdguess)
    [raval, decval] = solar_equ_coords(jdguess)
    lstval = LST(jdguess,longit)
    alt3 = getalt(raval, decval, lstval, lat)
    err = alt3 - alt
    if abs(err) < alterr:
        return jdguess
    deriv = (alt3 - alt2) / delta
    while( abs(err) > alterr and i < maxstep):
        jdguess = jdguess - err/deriv
        #[raval, decval] = SunCoords(jdguess)
        [raval, decval] = solar_equ_coords(jdguess)
        lstval = LST(jdguess, longit)
        alt3 = getalt(raval, decval, lstval, lat)
        err = alt3 - alt
        i = i + 1
        if i == maxstep-1:
            sys.stderr.write("Error determining the time at which the sun is at a given altitude.\n")
            #exit(2)
    if(i >= maxstep-1):
        jdguess = -1000.
    return jdguess

# Convert from UT calendar date to JD
def cal2jd(year_, month_, day_, hour_, minute_, second_):
    year = year_
    month = month_
    day = day_
    hour = hour_
    minute = minute_
    second = second_
    if(month <= 2):
        y = year - 1
        m = month + 12
    else:
        y = year
        m = month
    A = int(y / 100.)
    B = 2 - A + int(A/4.)
    jdout = int(365.25 * y) + int(30.6001*(m+1)) + day + 1720994.5
    jdout = jdout + (hour/24.) + (minute/1440.) + (second/86400.)
    if y > 1583:
        jdout = jdout + B
    return jdout

# Convert from JD to calendar date at specified time-zone
# tz = hours to add to UT time to get Local Time
def jd2cal(jd, tz):
    jdstore = jd + tz/24.0
    jdin = jdstore + 0.5
    z = math.floor(jdin)
    if z > 2299161:
        w = math.floor((4*z-7468865)/146097)
        a = 1+z+w-math.floor(w/4)
    else:
        a = z
    b = a+1524
    c = math.floor((20*b-2442)/7305)
    d = math.floor((1461*c)/4)
    e = math.floor(((b-d)*10000)/306001)
    da=int(b-d-z-math.floor(e*153/5)+jdin)
    if e <= 13:
        mo = int(e-1)
    else:
        mo = int(e-13)
    ye = int(math.floor(c-4716))
    if (mo <= 2):
        ye += 1
    ut = 24.0*(jdstore + 0.5 - math.floor(jdstore + 0.5))
    h = int(ut)
    m = int((ut - h) * 60.0)
    s = ((ut - h) * 60.0 - m)*60
    return [ye, mo, da, h, m, s]
    

# Determine the alt-degree twilight times (in JD) for a given night.
# nightyr, nightmo, and nightday are integers giving the year, month
# and day, lat is in degrees, longit is in degrees East.
# alt is in degrees
def gettwilightJD(alt, nightyr, nightmo, nightday, lat, longit):
    # First determine jd at approximate true midnight (not corrected
    # for equation of time) for this location.
    jdut = julian_day(nightyr, nightmo, nightday, 24, 0, 0)
    jdmidapprox = jdut - longit/15./24.
    
    # Check if the sun is always above or below the altitude threshold on
    # a given night.
    #[raval, decval] = SunCoords(jdmidapprox)
    [raval, decval] = solar_equ_coords(jdmidapprox)
    
    if lat >= 0.0:
        if 90.0 - lat + decval < alt:
            # The sun is never up, set the start and stop times to be
            # midnight +- 0.5
            return [jdmidapprox -0.5, jdmidapprox + 0.5]
        elif lat - (90.0 - decval) > alt:
            # The sun is always up, set the start and stop times to be
            # midnight
            return [jdmidapprox, jdmidapprox]
    else:
        if 90.0 + lat - decval < alt:
            return [jdmidapprox -0.5, jdmidapprox + 0.5]
        elif -lat - (90.0 + decval) > alt:
            # The sun is always up, set the start and stop times to be
            # midnight
            return [jdmidapprox, jdmidapprox]
        
    # Select trial evening and morning JD values
    ha = alt_to_HA(alt, decval, lat)
    jdeveapprox = jdmidapprox - ha/360.0
    jdeve = jd_sun_alt(alt, jdeveapprox, lat, longit)
    
    jdmornapprox = jdmidapprox + ha/360.0
    jdmorn = jd_sun_alt(alt, jdmornapprox, lat, longit)
    
    return [jdeve, jdmorn]

P_alt =[] #GET RID OF THESE THROUGH INCORPORATING THEM INTO THE TRANSITEVENT OBJECT
P_mag = []
P_flag = []
P_objprio = []
P_period = []
P_dur = []
Prios = []

periods = []
mags = []
durs = []

shallowmagmin = 14.948
shallowdurmin = 0.0128042
shallowdurmax = 0.694507
ptotal = 0.0

#  assign each object a priority based on magnitude, orbital period,
#  visibility flag, altitude at mid-transit, and number of subsequent
#  Bieryla nights
def pre_priority(flag, period, mag, altcent, objprio, HalfDuration):
    P_alt.append(abs(math.cos(math.radians(90 - altcent)))) # is the absolute value sign okay? 

    global ptotal

    mags.append(mag)

    if flag == 'OIBEO':
        P_flag.append(1.0)
    elif flag == '-IBEO' or flag == 'OIBE-':
        P_flag.append(0.75)
    elif flag == 'OIB--' or flag == '--BEO':
        P_flag.append(0.5)
    elif flag == 'OI---' or flag == '-IBE-' or flag == '---EO':
        P_flag.append(0.1)
    else:
        P_flag.append(0.0)

    P_objprio.append(10.0/(2.0**int(objprio)))

    ptotal += period
    periods.append(period)

    durs.append(HalfDuration*2.0)

# magnitude, period, and duration require comparison to all their values
# for proper priority assignment, so this is an addendum to account for that;
# additionally, it calculates relative priorities for each component listed
# def priority():
#     global shallowmagmin    # Needed to modify global copy of shallowmagmin
#     global shallowdurmin
#     global shallowdurmax

#     for i in range(0, len(mags)):
#         P_mag.append(10.0**(-0.02*(mags[i] - shallowmagmin))) 

#         P_period.append(periods[i]/(ptotal/len(periods))) 

#         P_dur.append(0.25 + ((0.5*(durs[i] - shallowdurmin))/(shallowdurmax - shallowdurmin))) 

#     minmag = (10.0**(-0.02*(shallowmagmin - shallowmagmin)))
#     maxmag = (10.0**(-0.02*(7.438 - shallowmagmin)))

#     minper = 0.114018/(ptotal/len(periods))
#     maxper = 32.9719/(ptotal/len(periods))

#     minalt = math.cos(math.radians(90 - 30.0))
#     maxalt = math.cos(math.radians(90 - 90.0))

#     mindur = 0.25 + ((0.5*(shallowdurmin - shallowdurmin))/(shallowdurmax - shallowdurmin))
#     maxdur = 0.25 + ((0.5*(shallowdurmax - shallowdurmin))/(shallowdurmax - shallowdurmin))

#     c_mag = math.log(1.5/0.5)/(math.log(maxmag/minmag))
#     c_alt = math.log(1.25/0.75)/(math.log(maxalt/minalt))
#     c_per = math.log(1.5/0.5)/(math.log(maxper/minper))
#     c_dur = math.log(1.1/0.9)/(math.log(maxdur/mindur))

#     for i in range(0, len(mags)):
#         P_mag[i] = 0.5 * ((P_mag[i]/minmag)**c_mag)
#         P_alt[i] = 0.75 * ((P_alt[i]/minalt)**c_alt)
#         P_period[i] = 0.5 * ((P_period[i]/minper)**c_per)
#         P_dur[i] = 0.9 * ((P_dur[i]/mindur)**c_dur)

#     for i in range(0, len(P_mag)):
#         Prios.append(P_alt[i] * P_mag[i] * P_flag[i] * P_objprio[i] * P_period[i] * P_dur[i])

# year, month, day are the dates to start and stop the check on
# RA and Dec are in degrees, period and HalfDuration are in days.
# Epoch is a JulianDate.
def FindVisibleTransits(yearstart, monthstart, daystart, yearstop, monthstop, daystop, RA, Dec, Period, Epoch, HalfDuration, obslat, obslongit, obstz, twialt, objalt, objname, objprio, objfutype, magv, starobj):
    # First generate a list of all transit start, stop, and midpoint times
    # Within the given date range.
    [jd1, jd2] = gettwilightJD(twialt, yearstart, monthstart, daystart, obslat, obslongit)
    if yearstart != yearstop or monthstart != monthstop or daystart != daystop:
        [jdtmp, jd2] = gettwilightJD(twialt, yearstop, monthstop, daystop, obslat, obslongit)
    nstart1 = math.floor((jd1 - Epoch + HalfDuration)/Period)
    nstop1 = math.floor((jd1 - Epoch - HalfDuration)/Period)
    if nstart1 > nstop1:
        # The object is in transit at evening twilight on the first night
        nstart = nstart1
    else:
        nstart = math.floor((jd1 - Epoch)/Period) + 1
    nstart2 = math.floor((jd2 - Epoch + HalfDuration)/Period)
    nstop2 = math.floor((jd2 - Epoch - HalfDuration)/Period)
    if nstart2 > nstop1:
        # The object is in transit at morning twilight on the last night
        nstop = nstart2
    else:
        nstop = math.floor((jd1 - Epoch)/Period)
    jdcent = Epoch + Period*numpy.array(range(int(nstart),int(nstop+1)))
    if len(jdcent) == 0:
        # No transits are in this range
        return []
    # Covert the series of transit center BJDs from the ephemeris into JDs
    command="vartools -i - -converttime input bjd output jd radec fix "+str(RA)+" "+str(Dec)+" -o - -quiet"
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
    s=""
    for jd in jdcent:
        s=s+"%.17g 0 0.1\n" % jd
    (output, dum) = process.communicate(s)
    jdcent_corr=[]
    for f in output.split("\n"):
        if len(f.split()) > 0:
            jdcent_corr.append(float(f.split()[0]))
    jdcent = numpy.array(jdcent_corr)
    transits = []
    # Now for each JD check if the object is above the horizon
    # and the sun is below the threshold for either the ingress,
    # midpoint, or egress
    ptotal = 0.0
    count = 0.0
    for jd in jdcent:
        jdOOTstart = jd - HalfDuration - 0.0208333
        jdstart = jd - HalfDuration
        jdstop = jd + HalfDuration
        jdOOTstop = jd + HalfDuration + 0.0208333
        [rasun, decsun] = solar_equ_coords(jdstart)
        lststart = LST(jdstart,obslongit)
        altsunstart = getalt(rasun, decsun, lststart, obslat)
        altstart = getalt(RA, Dec, lststart, obslat)
        [rasun, decsun] = solar_equ_coords(jd)
        lstcent = LST(jd,obslongit)
        altsuncent = getalt(rasun, decsun, lstcent, obslat)
        altcent = getalt(RA, Dec, lstcent, obslat)
        [rasun, decsun] = solar_equ_coords(jdstart)
        lststop = LST(jdstop,obslongit)
        altsunstop = getalt(rasun, decsun, lststop, obslat)
        altstop = getalt(RA, Dec, lststop, obslat)
        if altsunstart < twialt and altstart > objalt:
            flag = "I"
        else:
            flag = "-"
        if altsuncent < twialt and altcent > objalt:
            flag = flag+"B"
        else:
            flag = flag+"-"
        if altsunstop < twialt and altstop > objalt:
            flag = flag+"E"
        else:
            flag = flag+"-"
            
        if flag != "---":
            [rasun, decsun] = solar_equ_coords(jdOOTstart)
            lstOOTstart = LST(jdOOTstart,obslongit)
            altsunOOTstart = getalt(rasun, decsun, lstOOTstart, obslat)
            altOOTstart = getalt(RA, Dec, lstOOTstart, obslat)
            [rasun, decsun] = solar_equ_coords(jdOOTstop)
            lstOOTstop = LST(jdOOTstop,obslongit)
            altsunOOTstop = getalt(rasun, decsun, lstOOTstop, obslat)
            altOOTstop = getalt(RA, Dec, lstOOTstop, obslat)
            if altsunOOTstart < twialt and altOOTstart > objalt:
                flag = "O"+flag
            else:
                flag = "-"+flag
            if altsunOOTstop < twialt and altOOTstop > objalt:
                flag = flag+"O"
            else:
                flag = flag+"-"
            [yr, mo, da, hr, mi, se] = jd2cal(jd, obstz)
            if hr >= 12:
                yr1 = yr
                mo1 = mo
                da1 = da
                hr1 = hr
                mi1 = mi
                se1 = se
                [yr2, mo2, da2, hr2, mi2, se2] = jd2cal(jd+1.0, obstz)
            else:
                yr2 = yr
                mo2 = mo
                da2 = da
                hr2 = hr
                mi2 = mi
                se2 = se
                [yr1, mo1, da1, hr1, mi1, se1] = jd2cal(jd-1.0, obstz)
            night = "%04d.%02d.%02d/%02d" % (yr1, mo1, da1, da2)
            [yrsta, mosta, dasta, hrsta, mista, sesta] = jd2cal(jdstart, obstz)
            [yrsto, mosto, dasto, hrsto, misto, sesto] = jd2cal(jdstop, obstz)
            [ramoon, decmoon] = lunar_equ_coords(jd)
            moonsep = angular_separation(RA, Dec, ramoon, decmoon)
            ltstart = "%02d:%02d" % (hrsta, mista)
            ltcent = "%02d:%02d" % (hr, mi)
            ltstop = "%02d:%02d" % (hrsto, misto)
            if objfutype == "pri" or objfutype == "sec":
                # pre_priority(flag, Period, magv, altcent, objprio, HalfDuration)
                ptotal += Period
                count ++
                transits.append(TransitEvent(objname, objprio, math.floor((jd - Epoch + HalfDuration)/Period), night, flag, jdstart, ltstart, altstart, jd, ltcent, altcent, jdstop, ltstop, altstop, moonsep, objfutype, magv))
            elif objfutype == "odd" or objfutype == "even":
                objNtran = int(round((jd - Epoch)/Period)) % 2
                if objNtran == 0 and objfutype == "even":
                    priotoshow = objprio
                elif objNtran == 1 and objfutype == "odd":
                    priotoshow = objprio
                else:
                    priotoshow = "..."
                # pre_priority(flag, Period, magv, altcent, objprio, HalfDuration)
                ptotal += Period
                count ++
                transits.append(TransitEvent(objname, priotoshow, math.floor((jd - Epoch + HalfDuration)/Period), night, flag, jdstart, ltstart, altstart, jd, ltcent, altcent, jdstop, ltstop, altstop, moonsep, objfutype, magv))
    for obj in transits:
        obj.get_period(ptotal,count)
        obj.priority()
    return transits

# Run a shell command
def runcommand(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    os.waitpid(process.pid, 0)
    return process.stdout.read().strip()

# procedure to get the uname/pwd for transit.php
def getpasswdHATRED():
    # First check to see if the pwd file exists in the ~/.inspectcandidates
    # directory
    homedir=os.getenv("HOME")
    try:
        testfile = open(homedir+'/.inspectcandidates/passwd',"r")
    except IOError:
        # It doesn't exist, so prompt for the password
        root = Tk()
        dialog = userpwddialog(root,"HATRED")
        vals = dialog.activate()
        user = vals[0]
        pwd = vals[1]
        if not os.path.isdir(homedir+'/.inspectcandidates/'):
            runcommand('mkdir -p '+homedir+'/.inspectcandidates/')
        testfile2 = open(homedir+"/.inspectcandidates/passwd","w")
        testfile2.write(base64.b64encode(user) + "\n")
        testfile2.write(base64.b64encode(pwd) + "\n")
        testfile2.close()
        os.chmod(homedir+'/.inspectcandidates/passwd',(stat.S_IRUSR | stat.S_IWUSR))
    else:
        user = base64.b64decode(testfile.readline())
        pwd = base64.b64decode(testfile.readline())
        testfile.close()
    outvals=[user,pwd]
    return outvals

# procedure to get the uname/pwd for transit.php
def getpasswdHSCAND():
    # First check to see if the pwd file exists in the ~/.inspectcandidates
    # directory
    homedir=os.getenv("HOME")
    try:
        testfile = open(homedir+'/.inspectcandidates/hscandpasswd',"r")
    except IOError:
        # It doesn't exist, so prompt for the password
        root = Tk()
        dialog = userpwddialog(root,"HSCAND")
        vals = dialog.activate()
        user = vals[0]
        pwd = vals[1]
        if not os.path.isdir(homedir+'/.inspectcandidates/'):
            runcommand('mkdir -p '+homedir+'/.inspectcandidates/')
        testfile2 = open(homedir+"/.inspectcandidates/hscandpasswd","w")
        testfile2.write(base64.b64encode(user) + "\n")
        testfile2.write(base64.b64encode(pwd) + "\n")
        testfile2.close()
        os.chmod(homedir+'/.inspectcandidates/hscandpasswd',(stat.S_IRUSR | stat.S_IWUSR))
    else:
        user = base64.b64decode(testfile.readline())
        pwd = base64.b64decode(testfile.readline())
        testfile.close()
    outvals=[user,pwd]
    return outvals

def SecondaryParams(h, k, b2pri, Period, Tcpri, HalfDurpri):
    J = math.sqrt(1.0 - k*k - h*h)
    lpri = math.atan2((h+1.0-k*k/(1.0+J)),k+k*h/(1.0+J))-k*J/(1.0+h)
    lsec = math.atan2((h-1.0+k*k/(1.0+J)),k-k*h/(1+J))+k*J/(1.0-h)
    psec = (lsec - lpri)/(2.0*math.pi)
    while psec < 0: psec += 1.0
    while psec >= 1.0: psec -= 1.0
    esec = Tcpri + psec*Period
    b2sec = b2*((1.0+h)/(1.0-h))*((1.0+h)/(1.0-h))
    om = 1.0/(HalfDurpri)
    omsec = om*math.sqrt((1-b2)/(1-b2sec))*(1.0-h)/(1.0+h)
    HalfDursec = 1.0/omsec
    return [esec, HalfDursec]
    


# Look-up a candidate in the HATRED database
def QueryHATRED(obj):
    global isimportmysql
    global mode
    if isimportmysql == 0:
        import MySQLdb
        isimportmysql = 1
    [dbuser, dbpwd] = getpasswdHATRED()
    db = MySQLdb.connect(host='hat.astro.princeton.edu',user=dbuser,passwd=dbpwd,db='HATRED')
    db.cur = db.cursor()
    db.cur.execute('select HTRname, HTRra, HTRdec, HTRP, 2400000+HTRE+0.5*HTRP*HTRq, HTRq*0.5*HTRP, HTRTODO, HTRmagV from HTR where HTRname=+"'+obj+'"')
    if int(db.cur.rowcount) < 1:
        sys.stderr.write(obj+" is not in the HATRED database. Skipping\n")
        exit(3)
    row = db.cur.fetchone()
    prio="..."
    futype = "pri"
    if str(row[6]) != "None":
        for ftmp in row[6].split(','):
            if ftmp.split(':')[0] == "FLWO12" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                break
            elif ftmp.split(':')[0] == "FLWO12sec" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                futype = "sec"
                break
            elif ftmp.split(':')[0] == "FLWO12even" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                futype = "even"
                break
            elif ftmp.split(':')[0] == "FLWO12odd" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                futype = "odd"
                break
    if mode == "sec" or futype == "sec":
	# Check if this target has e/omega
	db.cur.execute('select HTRname, HTREP_h, HTREP_k, HTREP_B2 from HTREP where HTRname="'+obj+'" and HTREP_Status = 1 order by HTREP_PeakID asc')
	if int(db.cur.rowcount) > 0:
	    row2 = db.cur.fetchone()
	    if str(row2[1]) != "None" and str(row2[2]) != "None":
		h = float(row2[1])
		k = float(row2[2])
		if str(row2[3]) == "None":
		    b2 = 0.0
		else:
		    b2 = float(row2[3])
                [Tcsec, HalfDursec] = SecondaryParams(h, k, b2, float(row[3]), float(row[4]), float(row[5]))
                return ( row[1], row[2], row[3], Tcsec, HalfDursec, prio, futype, row[7] )
        # just increase the secondary time by 0.5*P
        Tcsec = float(row[4]) + 0.5*float(row[3])
        return ( row[1], row[2], row[3], Tcsec, row[5], prio, futype, row[7])
    vals = ( row[1], row[2], row[3], row[4], row[5], prio, futype, row[7] )
    return vals

# Look-up a candidate in the HATRED database
def QueryALLHATRED():
    global isimportmysql
    global mode
    if isimportmysql == 0:
        import MySQLdb
        isimportmysql = 1
    [dbuser, dbpwd] = getpasswdHATRED()
    db = MySQLdb.connect(host='hat.astro.princeton.edu',user=dbuser,passwd=dbpwd,db='HATRED')
    db.cur = db.cursor()
    db.cur.execute('select HTRname, HTRra, HTRdec, HTRP, 2400000+HTRE+0.5*HTRP*HTRq, HTRq*0.5*HTRP, HTRTODO, HTRmagV from HTR where HTRTODO is not NULL')
    objects = []
    rows = db.cur.fetchall()
    futype = "pri"
    for row in rows:
        prio="..."
        if str(row[3]) == "None" or str(row[4]) == "None" or str(row[5]) == "None":
            continue
        FUs = []
        for ftmp in row[6].split(','):
            if ftmp.split(':')[0] == "FLWO12" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                FUs.append(("pri",prio))
            elif ftmp.split(':')[0] == "FLWO12sec" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                FUs.append(("sec",prio))
            elif ftmp.split(':')[0] == "FLWO12even" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                FUs.append(("even",prio))
            elif ftmp.split(':')[0] == "FLWO12odd" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                FUs.append(("odd",prio))
        for (futype,prio) in FUs:
            if mode == "sec" or futype == "sec":
                # Check if this target has e/omega
                db.cur.execute('select HTRname, HTREP_h, HTREP_k, HTREP_B2 from HTREP where HTRname="'+obj+'" and HTREP_Status = 1 order by HTREP_PeakID asc')
                dumcheck = 0
                if int(db.cur.rowcount) > 0:
                    row2 = db.cur.fetchone()
                    if str(row2[1]) != "None" and str(row2[2]) != "None":
                        h = float(row2[1])
                        k = float(row2[2])
                        if str(row2[3]) == "None":
                            b2 = 0.0
                        else:
                            b2 = float(row2[3])
                        [Tcsec, HalfDursec] = SecondaryParams(h, k, b2, float(row[3]), float(row[4]), float(row[5]))
                        objects.append(Star(row[0], prio, float(row[1]), float(row[2]), float(row[3]), Tcsec, HalfDursec, futype, float(row[7])))
                        dumcheck=1
            # just increase the secondary time by 0.5*P
                if dumcheck == 0:
                    Tcsec = float(row[4]) + 0.5*float(row[3])
                    objects.append(Star(row[0], prio, float(row[1]), float(row[2]), float(row[3]), Tcsec, float(row[5]), futype, float(row[7])))
            else:
                objects.append(Star(row[0], prio, float(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]), futype, float(row[7])))
    return objects

# Look-up a candidate in the HSCAND database
def QueryHSCAND(obj):
    global isimportmysql
    global mode
    if isimportmysql == 0:
        import MySQLdb
        isimportmysql = 1
    [dbuser, dbpwd] = getpasswdHSCAND()
    db = MySQLdb.connect(host='hatsouth.astro.princeton.edu',user=dbuser,passwd=dbpwd,db='HSCAND')
    db.cur = db.cursor()
    db.cur.execute('select HATSname, HATSra, HATSdec, HATSP, 2400000+HATSE+0.5*HATSP*HATSq, HATSq*0.5*HATSP, HATSTODO, HATSmagV from HATS where HATSname=+"'+obj+'"')
    if int(db.cur.rowcount) < 1:
        sys.stderr.write(obj+" is not in the HSCAND database. Skipping\n")
        exit(3)
    row = db.cur.fetchone()
    prio="..."
    futype="pri"
    if str(row[6]) != "None":
	for ftmp in row[6].split(','):
            if ftmp.split(':')[0] == "PHFU" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                futype = "pri"
                break
            elif ftmp.split(':')[0] == "PHFUodd" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                futype = "odd"
                break
            elif ftmp.split(':')[0] == "PHFUeven" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                futype = "even"
                break
            elif ftmp.split(':')[0] == "PHFUsec" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                futype = "sec"
                break
    if mode == "sec" or futype == "sec":
	# Check if this target has e/omega
	db.cur.execute('select HATSname, HATSEP_h, HATSEP_k, HATSEP_B2 from HATSEP where HATSname="'+obj+'" and HATSEP_Status = 1 order by HATSEP_PeakID asc')
	if int(db.cur.rowcount) > 0:
	    row2 = db.cur.fetchone()
	    if str(row2[1]) != "None" and str(row2[2]) != "None":
		h = float(row2[1])
		k = float(row2[2])
		if str(row2[3]) == "None":
		    b2 = 0.0
		else:
		    b2 = float(row2[3])
                [Tcsec, HalfDursec] = SecondaryParams(h, k, b2, float(row[3]), float(row[4]), float(row[5]))
                return ( row[1], row[2], row[3], Tcsec, HalfDursec, prio, futype, row[7] )
        # just increase the secondary time by 0.5*P
        Tcsec = float(row[4]) + 0.5*float(row[3])
        return ( row[1], row[2], row[3], Tcsec, row[5], prio, futype, row[7])
    vals = ( row[1], row[2], row[3], row[4], row[5], prio, futype, row[7] )
    return vals

# Look-up a candidate in the HSCAND database
def QueryALLHSCAND():
    global isimportmysql
    global mode
    if isimportmysql == 0:
        import MySQLdb
        isimportmysql = 1
    [dbuser, dbpwd] = getpasswdHSCAND()
    db = MySQLdb.connect(host='hatsouth.astro.princeton.edu',user=dbuser,passwd=dbpwd,db='HSCAND')
    db.cur = db.cursor()
    db.cur.execute('select HATSname, HATSra, HATSdec, HATSP, 2400000+HATSE+0.5*HATSP*HATSq, HATSq*0.5*HATSP, HATSTODO, HATSmagV from HATS where HATSTODO is not NULL')
    objects = []
    rows = db.cur.fetchall()
    for row in rows:
        prio="..."
        futype="pri"
        FUs = []
        for ftmp in row[6].split(','):
            if ftmp.split(':')[0] == "PPHFU" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                FUs.append(("pri",prio))
            elif ftmp.split(':')[0] == "PPHFUodd" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                FUs.append(("odd",prio))
            elif ftmp.split(':')[0] == "PPHFUeven" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                FUs.append(("even",prio))
            elif ftmp.split(':')[0] == "PPHFUsec" and len(ftmp.split(':')) > 1:
                prio = ftmp.split(':')[1]
                FUs.append(("sec",prio))
            elif ftmp.split(':')[0] == "CPHFU" and len(ftmp.split(':')) > 1:
                prio = "C"+ftmp.split(':')[1]
                FUs.append(("pri",prio))
            elif ftmp.split(':')[0] == "CPHFUodd" and len(ftmp.split(':')) > 1:
                prio = "C"+ftmp.split(':')[1]
                FUs.append(("odd",prio))
            elif ftmp.split(':')[0] == "CPHFUeven" and len(ftmp.split(':')) > 1:
                prio = "C"+ftmp.split(':')[1]
                FUs.append(("even",prio))
            elif ftmp.split(':')[0] == "CPHFUsec" and len(ftmp.split(':')) > 1:
                prio = "C"+ftmp.split(':')[1]
                FUs.append(("sec",prio))
        if len(FUs) == 0:
            FUs.append(("pri","..."))
        for (futype,prio) in FUs:
            if mode == "sec" or futype == "sec":
                # Check if this target has e/omega
                db.cur.execute('select HATSname, HATSEP_h, HATSEP_k, HATSEP_B2 from HATSEP where HATSname="'+obj+'" and HATSEP_Status = 1 order by HATSEP_PeakID asc')
                dumcheck = 0
                if int(db.cur.rowcount) > 0:
                    row2 = db.cur.fetchone()
                    if str(row2[1]) != "None" and str(row2[2]) != "None":
                        h = float(row2[1])
                        k = float(row2[2])
                        if str(row2[3]) == "None":
                            b2 = 0.0
                        else:
                            b2 = float(row2[3])
                        [Tcsec, HalfDursec] = SecondaryParams(h, k, b2, float(row[3]), float(row[4]), float(row[5]))
                        objects.append(Star(row[0], prio, float(row[1]), float(row[2]), float(row[3]), Tcsec, HalfDursec, futype, float(row[7])))
                        dumcheck=1
                if dumcheck == 0:
                    # just increase the secondary time by 0.5*P
                    Tcsec = float(row[4]) + 0.5*float(row[3])
                    objects.append(Star(row[0], prio, float(row[1]), float(row[2]), float(row[3]), Tcsec, float(row[5]), futype, float(row[7])))
            else:
                objects.append(Star(row[0], prio, float(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]), futype, float(row[7])))
    return objects

# Look-up a candidate in the TEP database
def QueryTEP(obj):
    command = 'tepquery --id '+obj+' -p star_ra,star_dec,pl_period,pl_Tc,pl_T14,star_magV'
    val=runcommand(command)
    if val == "":
        exit(3)
    else:
        return val.split()

# Look-up a candidate in the TEP database
def QueryALLTEP():
    command = 'tepquery --all -p pl_name,star_ra,star_dec,pl_period,pl_Tc,pl_T14,star_magV'
    objects = []
    val=runcommand(command)
    if val == "":
        exit(3)
    else:
        for f in val.split('\n'):
            f2 = f.split()
            if str(f2[1]) == "None" or str(f2[2]) == "None" or str(f2[3]) == "None" or str(f2[4]) == "None" or str(f2[5]) == "None":
                continue
            if str(f2[6]) == "None":
                objects.append(Star(f2[0], "...", float(f2[1]), float(f2[2]), float(f2[3]), float(f2[4]), 0.5*float(f2[5]), "pri", 99.9))
            else:
                objects.append(Star(f2[0], "...", float(f2[1]), float(f2[2]), float(f2[3]), float(f2[4]), 0.5*float(f2[5]), "pri", float(f2[6])))
    return objects            

# Class to hold the dialog for prompting the user for the DB uname/pwd
class userpwddialog:
    def __init__(self, parent, DBname):
        self.dialog = Pmw.Dialog(parent,
                                 buttons = ('OK','Cancel'),
                                 title = DBname+' Database access',
                                 defaultbutton = 'Ok',
                                 buttonboxpos = 's',
                                 command = self.execute)
        self.dialog.withdraw()
        self.user = Pmw.EntryField(self.dialog.interior(),
                                   labelpos = 'w',
                                   value = '',
                                   label_text = 'username:')
        self.password = Pmw.EntryField(self.dialog.interior(),
                                       labelpos = 'w',
                                       value = '',
                                       label_text = 'password:',
                                       entry_show = '*')

        self.user.pack(side=TOP)
        self.password.pack(side=BOTTOM)
        self.passwordval = ""
        self.userval = ""
        

    def activate(self):
        self.dialog.activate()
        vals = [self.userval, self.passwordval]
        return vals

    def execute(self, result):
        if result == 'Cancel':
            exit(1)
        self.userval = self.user.getvalue()
        self.passwordval = self.password.getvalue()
        self.dialog.deactivate()

# procedure to get the hatred uname/pwd
def getdbuserpwd():
    # First check to see if the pwd file exists in the ~/.inspectcandidates 
    # directory
    homedir=os.getenv("HOME")
    try:
        testfile = open(homedir+'/.inspectcandidates/passwd',"r")
    except IOError:
        # It doesn't exist, so prompt for the password
        root = Tk()
        dialog = userpwddialog(root, "HATRED")
        vals = dialog.activate()
        dbuser = vals[0]
        dbpwd = vals[1]
        if not os.path.isdir(homedir+'/.inspectcandidates/'):
            runcommand('mkdir -p '+homedir+'/.inspectcandidates/')
        testfile2 = open(homedir+"/.inspectcandidates/passwd","w")
        testfile2.write(base64.b64encode(dbuser) + "\n")
        testfile2.write(base64.b64encode(dbpwd) + "\n")
        testfile2.close()
        os.chmod(homedir+'/.inspectcandidates/passwd',(stat.S_IRUSR | stat.S_IWUSR))
    else:
        dbuser = base64.b64decode(testfile.readline())
        dbpwd = base64.b64decode(testfile.readline())
        testfile.close()
    outvals=[dbuser,dbpwd]
    return outvals

# procedure to get the hscand uname/pwd
def gethscanddbuserpwd():
    # First check to see if the pwd file exists in the ~/.inspectcandidates 
    # directory
    homedir=os.getenv("HOME")
    try:
        testfile = open(homedir+'/.inspectcandidates/hscandpasswd',"r")
    except IOError:
        # It doesn't exist, so prompt for the password
        root = Tk()
        dialog = userpwddialog(root, "HSCAND")
        vals = dialog.activate()
        dbuser = vals[0]
        dbpwd = vals[1]
        if not os.path.isdir(homedir+'/.inspectcandidates/'):
            runcommand('mkdir -p '+homedir+'/.inspectcandidates/')
        testfile2 = open(homedir+"/.inspectcandidates/hscandpasswd","w")
        testfile2.write(base64.b64encode(dbuser) + "\n")
        testfile2.write(base64.b64encode(dbpwd) + "\n")
        testfile2.close()
        os.chmod(homedir+'/.inspectcandidates/hscandpasswd',(stat.S_IRUSR | stat.S_IWUSR))
    else:
        dbuser = base64.b64decode(testfile.readline())
        dbpwd = base64.b64decode(testfile.readline())
        testfile.close()
    outvals=[dbuser,dbpwd]
    return outvals

# procedure to start up db services
def startdb():
    vals = getdbuserpwd(parent)
    dbuser = vals[0]
    dbpwd = vals[1]
    db = MySQLdb.connect(host='hat',user=dbuser,passwd=dbpwd,db='HATRED')
    db.cur = db.cursor()
    return db

# procedure to start up hscanddb services
def starthscanddb():
    vals = gethscanddbuserpwd(parent)
    dbuser = vals[0]
    dbpwd = vals[1]
    hscanddb = MySQLdb.connect(host='hatsouth',user=dbuser,passwd=dbpwd,db='HSCAND')
    hscanddb.cur = hscanddb.cursor()
    return hscanddb

# parse the command line:
dateval=""
datestart=""
datestop=""
lon=""
lat=""
tz=""
site="flwo"
obj=""
name=""
ra=""
dec=""
per=""
Tc=""
q=""
mode="pri"
allHN=0
allHS=0
allTEP=0
TwiAlt=-12.0
ObjAlt=30.0
magv = "" 
try:
    opts, args = getopt.getopt(sys.argv[1:], "d:l:p:q:", ["start=", "stop=", "lon=", "lat=", "tz=", "obj=", "name=", "ra=", "dec=", "tc=", "secondary", "allHN", "allHS", "allTEP", "twialt=", "objalt=", "magv="])

except getopt.GetoptError:
    Usage()
    sys.exit(2)

for opt, arg in opts:
    if opt == "-d":
        dateval=arg
    elif opt == "--start":
        datestart=arg
    elif opt == "--stop":
        datestop=arg
    elif opt == "-l":
        site=arg
    elif opt == "--lon":
        lon=arg
    elif opt == "--lat":
        lat=arg
    elif opt == "--tz":
        tz=arg
    elif opt == "--obj":
        obj=arg
    elif opt == "--allHN":
        allHN=1
    elif opt == "--allHS":
        allHS=1
    elif opt == "--allTEP":
        allTEP=1
    elif opt == "--twialt":
        TwiAlt = arg
    elif opt == "--objalt":
        ObjAlt = arg
    elif opt == "--name":
        name=arg
    elif opt == "--ra":
        ra=arg
    elif opt == "--dec":
        dec=arg
    elif opt == "-p":
        per=arg
    elif opt == "--tc":
        Tc=arg
    elif opt == "-q":
        q=arg
    elif opt == "--magv":
        magv=arg
    elif opt == "--secondary":
        mode="sec"

if dateval == "" and (datestart == "" or datestop == ""):
    Usage()
    sys.exit(2)

if lon != "" and (lat == "" or tz == ""):
    Usage()
    sys.exit(2)

if lat != "" and (lon == "" or tz == ""):
    Usage()
    sys.exit(2)

if tz != "" and (lon == "" or lat == ""):
    Usage()
    sys.exit(2)

if lon == "" and lat == "" and tz == "":
    if site == "flwo":
        lon="-110.878"
        lat="31.681"
        tz="-7"
    elif site == "keck":
        lon="-155.477"
        lat="19.824"
        tz="-10"
    elif site == "ftn":
        lon="-156.2558"
        lat="20.7075"
        tz="-10"
    elif site == "fts":
	lon="149.0702778"
        lat="-31.2731944"
        tz="10"
    elif site == "lapalma":
        lon="-17.87916"
        lat="28.7583333"
        tz="0"
    elif site == "wise":
        lon="34.7633333"
        lat="30.5958333"
        tz="2"
    elif site == "piszkes":
        lon="19.895"
        lat="47.918"
        tz="2"
    elif site == "lasilla":
        lon="-70.729167"
        lat="-29.2566667"
        tz="-4"
    elif site == "ctio":
        lon="-70.816666499999997"
        lat="-30.1650000"
        tz="-4"
    elif site == "ohp":
        lon="5.7133333"
        lat="43.9308333"
        tz="0"
    elif site == "mcdonald":
        lon="-104.0225"
        lat="30.6714"
        tz="-6"
    elif site == "saao":
        lon="18.4776"
        lat="-33.9347"
        tz="-2"
    else:
        Usage()
        sys.exit(2)

if not (ra == "" and dec == "" and per == "" and Tc == "" and q == "" and magv == "") and \
        not (ra != "" and dec != "" and per != "" and Tc != "" and q != "" and magv != ""):
    Usage()
    sys.exit(2)

if dateval != "":
    yrstart=int(dateval.split('.')[0])
    mostart=int(dateval.split('.')[1][0:2])
    daystart=int(dateval.split('.')[1][2:])
    yrstop = int(yrstart)
    mostop = int(mostart)
    daystop = int(daystart)
else:
    yrstart=int(datestart.split('.')[0])
    mostart=int(datestart.split('.')[1][0:2])
    daystart=int(datestart.split('.')[1][2:])
    yrstop = int(datestop.split('.')[0])
    mostop = int(datestop.split('.')[1][0:2])
    daystop = int(datestop.split('.')[1][2:])

Objects = []
# Get objects from the Database if requested:
if name != "":
    HalfDuration = float(per)*0.5*float(q)
    Objects.append(Star(name, "...", float(ra), float(dec), float(per), float(Tc), HalfDuration, "pri", float(magv)))

if obj != "":
    if not re.search('^HTR',obj) and not re.search('^HATS[0-9][0-9][0-9]\-[0-9][0-9][0-9]$',obj):
        if mode == "sec":
            sys.stderr.write("Error: --secondary option not yet supported for TEPs")
            exit(3)
        vals = QueryTEP(obj)
        ra = vals[0]
        dec = vals[1]
        per = vals[2]
        Tc = vals[3]
        halfdur = float(vals[4])*0.5
        prio = "..."
        magv = vals[7]
        Objects.append(Star(obj, prio, float(ra), float(dec), float(per), float(Tc), float(halfdur), "pri", float(magv)))
    elif re.search('^HTR',obj):
        vals = QueryHATRED(obj)
        ra = vals[0]
        dec = vals[1]
        per = vals[2]
        Tc = vals[3]
        halfdur = float(vals[4])
        prio = vals[5]
        futype = vals[6]
        magv = vals[7]
        Objects.append(Star(obj, prio, float(ra), float(dec), float(per), float(Tc), float(halfdur), futype, float(magv)))
    elif re.search('^HATS[0-9][0-9][0-9]\-[0-9][0-9][0-9]$',obj):
        vals = QueryHSCAND(obj)
        ra = vals[0]
        dec = vals[1]
        per = vals[2]
        Tc = vals[3]
        halfdur = float(vals[4])
        prio = vals[5]
        futype = vals[6]
        magv = vals[7]
        Objects.append(Star(obj, prio, float(ra), float(dec), float(per), float(Tc), float(halfdur), futype, float(magv)))
    else:
        sys.stderr.write("Unknown object: "+obj+"\n")
        exit(2)
    
if allHN == 1:
    Objects.extend(QueryALLHATRED())

if allHS == 1:
    Objects.extend(QueryALLHSCAND())

if allTEP == 1:
    if mode == "sec":
        sys.stderr.write("Error: --secondary option not yet supported for TEPs")
        exit(3)
    Objects.extend(QueryALLTEP())

alltransits = []
for st in Objects:
    alltransits.extend(FindVisibleTransits(yrstart, mostart, daystart, yrstop, mostop, daystop, st.RA, st.Dec, st.P, st.Epoch, st.HalfDuration, float(lat), float(lon), float(tz), float(TwiAlt), float(ObjAlt), st.obj, st.prio, st.futype, st.magv))

# priority()
# for i in range(0, len(alltransits)):
#     alltransits[i].eventPriority = Prios[i]

alltransits_sort = sorted(alltransits, key=lambda trans: trans.jdcent)

for trans in alltransits_sort:
    trans.printvalues()
