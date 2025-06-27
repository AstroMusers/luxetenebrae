import astropy.units as u
from astropy import constants as const
import numpy as np
import logging
import sys, os

def is_ms_bh_pair(st1, st2):
    return ((st1 <= 9 and st2 == 14) or (st2 <= 9 and st1 == 14))

def is_bh_bh_pair(st1, st2):
    return (st1 == 14 and st2 == 14)

def check_merger_MT_hist(mt_history):
    return (mt_history[-1] == 6)

def check_merger_manual(r1, r2, sa, e):
    if e <= 1:
        return ((r1 + r2) >= sa)
    else:
        return False
    
def which_bh(stellar_type_1, stellar_type_2):
    if stellar_type_1 == 14:
        return 1
    elif stellar_type_2 == 14:
        return 0
    
def type_change_MT(preMT1, pstMT1, preMT2, pstMT2, start_time):
    maskchangeMT1 = np.zeros(len(preMT1), dtype='bool')
    maskchangeMT2 = np.zeros(len(preMT2), dtype='bool')

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
    
def extract_masked_systems_by_seed(Data, survivors):
    """
    Returns a dictionary mapping each seed to a dictionary of lists of all relevant entries.
    """
    SP = Data['BSE_System_Parameters']
    MT = Data['BSE_RLOF']
    CE = Data['BSE_Common_Envelopes']

    # print(SP.keys())
    seedsSP = SP['SEED'][()]
    seedsMT = MT['SEED'][()]
    seedsCE = CE['SEED'][()]

    seedsSP_mask = np.isin(seedsSP, survivors) 
    seedsMT_mask = np.isin(seedsMT, survivors)
    seedsCE_mask = np.isin(seedsCE, survivors)

    statusSP = SP['Evolution_Status'][()]

    stellarTypeZamsSP1   =  SP['Stellar_Type@ZAMS(1)'][()]
    stellarTypeZamsSP2   =  SP['Stellar_Type@ZAMS(2)'][()]
    stellarTypeSP1   =  SP['Stellar_Type(1)'][()]
    stellarTypeSP2   =  SP['Stellar_Type(2)'][()]  

    massZamsSP1 = SP['Mass@ZAMS(1)'][()] 
    massZamsSP2 = SP['Mass@ZAMS(2)'][()]

    semimajorAxisZamsSP = SP['SemiMajorAxis@ZAMS'][()] # in AU
    eccentricityZamsSP = SP['Eccentricity@ZAMS'][()]
    mergersSP = SP['Merger'][()]

    MT = Data['BSE_RLOF']

    seedsMT = MT['SEED'][()]

    eventsMT = MT['MT_Event_Counter'][()]

    stellarTypepreMT1   =  MT['Stellar_Type(1)<MT'][()]
    stellarTypepreMT2   =  MT['Stellar_Type(2)<MT'][()]  
    stellarTypepstMT1   =  MT['Stellar_Type(1)>MT'][()]
    stellarTypepstMT2   =  MT['Stellar_Type(2)>MT'][()] 

    masspreMT1 = MT['Mass(1)<MT'][()] 
    masspreMT2 = MT['Mass(2)<MT'][()] 
    masspstMT1 = MT['Mass(1)>MT'][()] 
    masspstMT2 = MT['Mass(2)>MT'][()]

    semimajorAxispreMT = MT['SemiMajorAxis<MT'][()] # in R_sun
    semimajorAxispstMT = MT['SemiMajorAxis>MT'][()]

    eccentricitypreMT = MT['Eccentricity<MT'][()]
    eccentricitypstMT = MT['Eccentricity>MT'][()]

    semimajorAxispreMT = (semimajorAxispreMT*const.R_sun).to(u.au).value
    semimajorAxispstMT = (semimajorAxispstMT*const.R_sun).to(u.au).value

    timepreMT = MT['Time<MT'][()]
    timepstMT = MT['Time>MT'][()]

    radiuspreMT1 = MT['Radius(1)<MT'][()]
    radiuspstMT1 = MT['Radius(1)>MT'][()] 
    radiuspreMT2 = MT['Radius(2)<MT'][()]
    radiuspstMT2 = MT['Radius(2)>MT'][()]

    CEafterMT = MT['CEE>MT'][()]

    CE = Data['BSE_Common_Envelopes']

    seedsCE = CE['SEED'][()]
    eventsCE = CE['CE_Event_Counter'][()]
    
    stellarTypepreCE1   =  CE['Stellar_Type(1)<CE'][()]
    stellarTypepreCE2   =  CE['Stellar_Type(2)<CE'][()]  
    stellarTypepstCE1   =  CE['Stellar_Type(1)'][()]
    stellarTypepstCE2   =  CE['Stellar_Type(2)'][()] 
    # PreCE refers to the mass just before the first CE event, while pstCE refers to the mass just after the last CE event.
    masspreCE1 = CE['Mass(1)<CE'][()] 
    masspreCE2 = CE['Mass(2)<CE'][()] 
    masspstCE1 = CE['Mass(1)>CE'][()] 
    masspstCE2 = CE['Mass(2)>CE'][()]

    semimajorAxispreCE = CE['SemiMajorAxis<CE'][()] # in R_sun
    semimajorAxispstCE = CE['SemiMajorAxis>CE'][()]

    semimajorAxispreCE = (semimajorAxispreCE*const.R_sun).to(u.au).value
    semimajorAxispstCE = (semimajorAxispstCE*const.R_sun).to(u.au).value

    timeCE = CE['Time'][()]
    circularizationTime = CE['Tau_Circ'][()]
    
    radiuspreCE1 = CE['Radius(1)<CE'][()] # in R_sun
    radiuspstCE1 = CE['Radius(1)>CE'][()] 
    radiuspreCE2 = CE['Radius(2)<CE'][()]
    radiuspstCE2 = CE['Radius(2)>CE'][()] 

    # Apply mask
    selected_seeds = seedsSP[seedsSP_mask]
    selected_status = statusSP[seedsSP_mask]
    selected_stellarTypeZamsSP1 = stellarTypeZamsSP1[seedsSP_mask]
    selected_stellarTypeZamsSP2 = stellarTypeZamsSP2[seedsSP_mask]
    selected_stellarTypeSP1 = stellarTypeSP1[seedsSP_mask]
    selected_stellarTypeSP2 = stellarTypeSP2[seedsSP_mask]
    selected_massZamsSP1 = massZamsSP1[seedsSP_mask]
    selected_massZamsSP2 = massZamsSP2[seedsSP_mask]
    selected_semimajorAxisZamsSP = semimajorAxisZamsSP[seedsSP_mask]  # in AU
    selected_eccentricityZamsSP = eccentricityZamsSP[seedsSP_mask]
    selected_mergersSP = mergersSP[seedsSP_mask]

    selected_seedsMT = seedsMT[seedsMT_mask]
    selected_eventsMT = eventsMT[seedsMT_mask]
    selected_stellarTypepreMT1 = stellarTypepreMT1[seedsMT_mask]
    selected_stellarTypepreMT2 = stellarTypepreMT2[seedsMT_mask]
    selected_stellarTypepstMT1 = stellarTypepstMT1[seedsMT_mask]
    selected_stellarTypepstMT2 = stellarTypepstMT2[seedsMT_mask]
    selected_masspreMT1 = masspreMT1[seedsMT_mask]
    selected_masspreMT2 = masspreMT2[seedsMT_mask]
    selected_masspstMT1 = masspstMT1[seedsMT_mask]
    selected_masspstMT2 = masspstMT2[seedsMT_mask]
    selected_semimajorAxispreMT = semimajorAxispreMT[seedsMT_mask]  # in AU
    selected_semimajorAxispstMT = semimajorAxispstMT[seedsMT_mask]  # in AU
    selected_eccentricitypreMT = eccentricitypreMT[seedsMT_mask]
    selected_eccentricitypstMT = eccentricitypstMT[seedsMT_mask]
    selected_timepreMT = timepreMT[seedsMT_mask]
    selected_timepstMT = timepstMT[seedsMT_mask]
    selected_radiuspreMT1 = radiuspreMT1[seedsMT_mask]  # in R_sun
    selected_radiuspstMT1 = radiuspstMT1[seedsMT_mask]  # in R_sun
    selected_radiuspreMT2 = radiuspreMT2[seedsMT_mask]  # in R_sun
    selected_radiuspstMT2 = radiuspstMT2[seedsMT_mask]  # in R_sun
    selected_ceafterMT = CEafterMT[seedsMT_mask]

    selected_seedsCE = seedsCE[seedsCE_mask]
    selected_eventsCE = eventsCE[seedsCE_mask]
    selected_stellarTypepreCE1 = stellarTypepreCE1[seedsCE_mask]
    selected_stellarTypepreCE2 = stellarTypepreCE2[seedsCE_mask]
    selected_stellarTypepstCE1 = stellarTypepstCE1[seedsCE_mask]
    selected_stellarTypepstCE2 = stellarTypepstCE2[seedsCE_mask]
    selected_masspreCE1 = masspreCE1[seedsCE_mask]
    selected_masspreCE2 = masspreCE2[seedsCE_mask]
    selected_masspstCE1 = masspstCE1[seedsCE_mask]
    selected_masspstCE2 = masspstCE2[seedsCE_mask]
    selected_semimajorAxispreCE = semimajorAxispreCE[seedsCE_mask]  # in AU
    selected_semimajorAxispstCE = semimajorAxispstCE[seedsCE_mask]  # in AU
    selected_timeCE = timeCE[seedsCE_mask]
    selected_circularizationTime = circularizationTime[seedsCE_mask]
    selected_radiuspreCE1 = radiuspreCE1[seedsCE_mask]  # in R_sun
    selected_radiuspstCE1 = radiuspstCE1[seedsCE_mask]  # in R_sun
    selected_radiuspreCE2 = radiuspreCE2[seedsCE_mask]  # in R_sun
    selected_radiuspstCE2 = radiuspstCE2[seedsCE_mask]  # in R_sun
    # Create a dictionary to hold the data for each seed
    # Build dictionary: seed -> dict of lists
    systems_by_seed_SP = {}
    systems_by_seed_MT = {}
    systems_by_seed_CE = {}

    for i in range(len(selected_seeds)):
        seed = selected_seeds[i]
        if seed not in systems_by_seed_SP:
            systems_by_seed_SP[seed] = {
                'SEED': seed,
                'Evolution_Status': selected_status[i],
                'Stellar_Type@ZAMS(1)': selected_stellarTypeZamsSP1[i],
                'Stellar_Type@ZAMS(2)': selected_stellarTypeZamsSP2[i],
                'Stellar_Type(1)': selected_stellarTypeSP1[i],
                'Stellar_Type(2)': selected_stellarTypeSP2[i],
                'Mass@ZAMS(1)': selected_massZamsSP1[i],
                'Mass@ZAMS(2)': selected_massZamsSP2[i],
                'SemiMajorAxis@ZAMS': selected_semimajorAxisZamsSP[i],
                'Eccentricity@ZAMS': selected_eccentricityZamsSP[i],
                'Merger': selected_mergersSP[i],
            }

    for i in range(len(selected_seedsMT)):
        seed = selected_seedsMT[i]
        if seed not in systems_by_seed_MT:
            systems_by_seed_MT[seed] = {
                'SEED_MT': seedsMT[i],
                'MT_Event_Counter': selected_eventsMT[i],
                'Stellar_Type(1)<MT': selected_stellarTypepreMT1[i],
                'Stellar_Type(2)<MT': selected_stellarTypepreMT2[i],
                'Stellar_Type(1)>MT': selected_stellarTypepstMT1[i],
                'Stellar_Type(2)>MT': selected_stellarTypepstMT2[i],
                'Mass(1)<MT': selected_masspreMT1[i],
                'Mass(2)<MT': selected_masspreMT2[i],
                'Mass(1)>MT': selected_masspstMT1[i],
                'Mass(2)>MT': selected_masspstMT2[i],
                'SemiMajorAxis<MT': selected_semimajorAxispreMT[i],  #
                'SemiMajorAxis>MT': selected_semimajorAxispstMT[i],  # in AU
                'Eccentricity<MT': selected_eccentricitypreMT[i],
                'Eccentricity>MT': selected_eccentricitypstMT[i],
                'Time<MT': selected_timepreMT[i],
                'Time>MT': selected_timepstMT[i],
                'CE_History': selected_ceafterMT[i],  # CE history
                'Radius(1)<MT': selected_radiuspreMT1[i],  # in R
                'Radius(1)>MT': selected_radiuspstMT1[i],  # in R_sun
                'Radius(2)<MT': selected_radiuspreMT2[i],  # in R
                'Radius(2)>MT': selected_radiuspstMT2[i],  # in R_sun
            }

    for i in range(len(selected_seedsCE)):
        seed = selected_seedsCE[i]
        if seed not in systems_by_seed_CE:
            systems_by_seed_CE[seed] = {
                'SEED_CE': seedsCE[i],
                'CE_Event_Counter': selected_eventsCE[i],
                'Stellar_Type(1)<CE': selected_stellarTypepreCE1[i],
                'Stellar_Type(2)<CE': selected_stellarTypepreCE2[i],
                'Stellar_Type(1)': selected_stellarTypepstCE1[i],  # Post CE
                'Stellar_Type(2)': selected_stellarTypepstCE2[i],  # Post CE
                'Mass(1)<CE': selected_masspreCE1[i],  # Pre CE
                'Mass(2)<CE': selected_masspreCE2[i],  # Pre CE
                'Mass(1)>CE': selected_masspstCE1[i],  # Post CE
                'Mass(2)>CE': selected_masspstCE2[i],  # Post CE
                'SemiMajorAxis<CE': selected_semimajorAxispreCE[i],  # in AU
                'SemiMajorAxis>CE': selected_semimajorAxispstCE[i],  # in AU
                'Time_CE': selected_timeCE[i],  # Time of the CE event
                "Tau_Circ": selected_circularizationTime[i],  # Circularization time
                'Radius(1)<CE': selected_radiuspreCE1[i],  # in R_sun
                'Radius(1)>CE': selected_radiuspstCE1[i],  # in R_sun
                'Radius(2)<CE': selected_radiuspreCE2[i],  # in R_sun
                'Radius(2)>CE': selected_radiuspstCE2[i],  # in R_sun
            }

    return systems_by_seed_SP, systems_by_seed_MT, systems_by_seed_CE

class StreamToLogger:
    def __init__(self, logger, log_level):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())

    def flush(self):
        pass