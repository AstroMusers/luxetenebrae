import numpy as np
from datetime import datetime
import os
import logging
import datetime as dt
import sys
from astropy import constants as const
from astropy import units as u
from scipy.integrate import quad
from scipy.interpolate import interp1d

pathToData = '/data/a.saricaoglu/repo/COMPAS'

def orbital_period(m1, m2, semimajax):
    # G =  39.4769264 # gravitational constant in AU^3 / (year^2 x Msun) 
    # M = (np.asarray(m1) + np.asarray(m2))
    # A = np.asarray(semimajax)
    # orbitalPeriod = (np.sqrt((4 * np.pi**2 * A**3)) / np.sqrt(G * M)) * 365.25
    G = const.G # units of m^3 kg^-1 s^-2
    m1 = m1 * const.M_sun # units of kg
    m1 = m1.to(u.kg)
    m2 = m2 * const.M_sun # units of kg
    m2 = m2.to(u.kg)
    A = semimajax * u.au # units of m
    A = A.to(u.m)
    orbitalPeriod = (np.sqrt((4 * np.pi**2 * A**3)) / np.sqrt(G * (m1 + m2))) # units of s
    orbitalPeriod = orbitalPeriod.to(u.hour).value
    print('max orbital period in hours:', np.max(orbitalPeriod))
    orbitalPeriod = orbitalPeriod / 24 # units of days
    print('max orbital period in days:', np.max(orbitalPeriod))

    return orbitalPeriod

def orbital_inclination(r1,r2,semimajax):

    #both radius and semimajax are in solar radii, no need to convert
    R = np.maximum(np.asarray(r1),np.asarray(r2))
    #print(R)
    RoverA = R/np.asarray(semimajax)

    return RoverA

def periapsis(semimajax, eccentricity):
    return (semimajax * (1 - eccentricity))

def searchability(orbital_inclination):

    inc = np.random.uniform(-np.pi,np.pi,len(orbital_inclination))
    cosi = np.absolute(np.cos(inc))
    print(cosi)
    searchability_index = [cosi <= np.cos(np.deg2rad(89.9))][0]

    return searchability_index

def envelope_ejection(r1,r2,postcesemimajax):
    R = np.asarray(r1) + np.asarray(r2)
    ejection = [R < np.asarray(postcesemimajax)][0]

    return ejection

def type_change_MT(preMT1, pstMT1, preMT2, pstMT2, start_time):
    maskchangeMT1 = np.zeros(len(preMT1), dtype='bool')
    maskchangeMT2 = np.zeros(len(preMT2), dtype='bool')
    s = datetime.now()

    #print(f"Run : {start_time}")

    for i in range(0, len(preMT1)):
        if (preMT1[i] == pstMT1[i]):
            maskchangeMT1[i] = False
            ##print(f'StellarType(1)<MT : {preMT1[i]}, StellarType(1)>MT : {pstMT1[i]} --> No change in stellar type during mass transfer event, maskchangeMT1[ {i} = {maskchangeMT1[i]}')

        else:
            maskchangeMT1[i] = True
            #print(f'StellarType(1)<MT : {preMT1[i]}, StellarType(1)>MT : {pstMT1[i]} --> Change in stellar type during mass transfer event, maskchangeMT1[ {i} = {maskchangeMT1[i]}')

    #print(f"Run : {start_time}")

    for i in range(0, len(preMT2)):
        if (preMT2[i] == pstMT2[i]):
            maskchangeMT2[i] = False
            ##print(f'StellarType(2)<MT : {preMT2[i]}, StellarType(2)>MT : {pstMT2[i]} --> No change in stellar type during mass transfer event, maskchangeMT1[ {i} = {maskchangeMT2[i]}')

        else:
            maskchangeMT2[i] = True
            #print(f'StellarType(2)<MT : {preMT2[i]}, StellarType(2)>MT : {pstMT2[i]} --> Change in stellar type during mass transfer event, maskchangeMT1[ {i} = {maskchangeMT2[i]}')


    return maskchangeMT1, maskchangeMT2

def type_change_CE(preCE1, pstCE1, preCE2, pstCE2, start_time):
    maskchangeCE1 = np.zeros(len(preCE1), dtype='bool')
    maskchangeCE2 = np.zeros(len(preCE2), dtype='bool')

    for i in range(0, len(preCE1)):
        if (preCE1[i] == pstCE1[i]):
            maskchangeCE1[i] = False
            ##print(f'StellarType(1)<CE : {preCE1[i]}, StellarType(1)>CE : {pstCE1[i]} --> No change in stellar type during mass transfer event, maskchangeCE1[ {i} = {maskchangeCE1[i]}')

        else:
            maskchangeCE1[i] = True
            #print(f'StellarType(1)<CE : {preCE1[i]}, StellarType(1)>CE : {pstCE1[i]} --> Change in stellar type during mass transfer event, maskchangeCE1[ {i} = {maskchangeCE1[i]}')

    #print(f"Run : {start_time}")

    for i in range(0, len(preCE2)):
        if (preCE2[i] == pstCE2[i]):
            maskchangeCE2[i] = False
            ##print(f'StellarType(2)<CE : {preCE2[i]}, StellarType(2)>CE : {pstCE2[i]}, --> No change in stellar type during mass transfer event, maskchangeCE1[ {i} = {maskchangeCE2[i]}')

        else:
            maskchangeCE2[i] = True
            #print(f'StellarType(2)<CE : {preCE2[i]}, StellarType(2)>CE : {pstCE2[i]} --> Change in stellar type during mass transfer event, maskchangeCE1[ {i} = {maskchangeCE2[i]}')

    return maskchangeCE1, maskchangeCE2


def find_last_mt(MTs):

    masklastMT = np.zeros(len(MTs), dtype='bool')
    last = 0
    for i in range(0, len(MTs)):
        if (MTs[i] == last):
            masklastMT[i-1] = False
            masklastMT[i] = True
        else:
            masklastMT[i] = True
            last = MTs[i]

    return masklastMT

def density(counts, bins):
    """
    Calculate the density of points in a histogram.
    
    Parameters:
    counts (array-like): The counts of points in each bin.
    bins (array-like): The edges of the bins.
    
    Returns:
    array: The density of points in each bin.
    """
    bin_widths = np.diff(bins)
    return counts / bin_widths

def kroupa_imf_continuous(m):
    """
    Kroupa IMF for continuous mass distribution.
    
    Parameters:
    m (float): Mass in solar masses.
    
    Returns:
    float: IMF value for the given mass.
    """
    if m <= 0.08:
        return m**(-0.3)
    elif m < 0.5:
        return m**(-1.3)
    else:
        return m**(-2.3)

def kroupa_imf(masses):
    """
    Calculate the Kroupa IMF for a given array of masses.
    
    Parameters:
    masses (array-like): Array of masses in solar masses.
    
    Returns:
    array: Kroupa IMF values for the given masses.
    """
    return np.array([kroupa_imf_continuous(m) for m in masses])

def kroupa_imf_normalized(m_min, m_max):
    """
    Normalize the Kroupa IMF for a given array of masses.
    
    Parameters:
    masses (array-like): Array of masses in solar masses.

    Returns:
    array: Normalized Kroupa IMF values for the given masses.
    """

    masses = np.linspace(m_min, m_max, 1000)
    normalization, _ = quad(kroupa_imf_continuous, m_min, m_max)
    imf = kroupa_imf(masses)
    return imf / normalization

# s = searchability(np.random.uniform(-1,1,10))

#print(np.shape(s))
#print(s)
print(const.G)
print(const.M_sun)
print(u.M_sun)
print(const.M_sun.to(u.kg))
print(u.M_sun.to(u.kg))
print((93*const.R_sun).to(u.au))
G = const.G # units of m^3 kg^-1 s^-2
m1 = (3*u.M_sun).to(u.kg) # units of kg
m2 = (8*u.M_sun).to(u.kg)  # units of kg
A = (7 * u.au).to(u.m) # units of m
orbitalPeriod = (np.sqrt((4 * np.pi**2 * A**3)) / np.sqrt(G * (m1 + m2))) 
print(orbitalPeriod)
orbinc = orbital_inclination([1,2,3],[4,5,6],[70,80,90])
print(orbinc)
se = searchability(np.random.uniform(-1,1,1000000))
print(np.sum(se))
print('percentage of searchability:', np.sum(se)/len(se))
print(np.arccos(2.228e-6))
print(np.deg2rad(89.99))
print(np.cos(np.deg2rad(89.99)))