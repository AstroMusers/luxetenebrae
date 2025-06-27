#### PARENT SCRIPT: secundus_processing_020525
#----> this version, instead of only MT data frame, reads CE too.
import os, sys
import numpy as np               # for handling arrays
import numpy.ma as ma
import pandas as pd
import logging            # for reading the COMPAS data
import time as t                      # for finding computation time
import datetime as dt
import h5py as h5                # for reading the COMPAS data
import matplotlib.pyplot as plt  # for plotting
import matplotlib
from luxetenebrae import utils as utils      # functions from utils.py
import astropy.units as u
from astropy import constants as const
# Import COMPAS specific scripts
import astropy.io.fits as fits
from astropy.table import Table

now = dt.datetime.now()
# Choose the mode to process
mode = 'WD_Enabled_Detailed'

# Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print("COMPAS_ROOT_DIR:", compasRootDir)
print(sys.path)

pathToData = os.path.join(compasRootDir, "luxetenebrae/runs", mode)
pathToFiles = os.path.join(compasRootDir, "luxetenebrae/files", mode)
pathToLogs = os.path.join(compasRootDir, "luxetenebrae/logs",mode)
# Create output directory if it doesn't exist
os.makedirs(pathToFiles, exist_ok=True)
os.makedirs(pathToLogs, exist_ok=True)

# Get the script name
script_name = os.path.basename(__file__)
# Configure logging
log_filename = f"{pathToLogs}/{script_name}_{now.strftime('%m.%d')}_{now.strftime('%H%M')}.log"
os.makedirs(os.path.dirname(log_filename), exist_ok=True)
logging.basicConfig(filename=log_filename, level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Redirect stdout and stderr to the log file
sys.stdout = utils.StreamToLogger(logging.getLogger('STDOUT'), logging.INFO)
sys.stderr = utils.StreamToLogger(logging.getLogger('STDERR'), logging.ERROR)
# Log the start of the script with the script name
logging.info(f'Script {script_name} started')

# Displays Time
s = dt.datetime.now()
starttime = t.ctime(t.time())
start = t.process_time()
start_time = s.strftime("%d%m%y") + "_" + s.strftime('%H%M')
print("Start time :", start_time)


# Choose an output hdf5 file to work with

ps = open(f'{compasRootDir}/luxetenebrae/files/survivors.txt', 'r')
survivors = ps.read().splitlines()


runs= [x[0] for x in os.walk(pathToData) if "Detailed_Output" in x[1]]
print('runs:', runs)
data_outputs_1 = []
data_outputs_2 = []
for run in runs:  # Adjust the slice as needed
    out1 = [f for f in os.listdir(run + '/Detailed_Output') if ".h5" in f]
    out2 = [f for f in os.listdir(run) if ".h5" in f]
    for f in out1:
        print("Reading file: ", run + '/Detailed_Output' + "/" + f)
        try:
            data = h5.File(run + '/Detailed_Output' + "/" + f)
            data_outputs_1.append(data)
        except:
            continue
    for f in out2:
        print("Reading file: ", run  + "/" + f)
        data = h5.File( run + "/" + f, 'r')
        data_outputs_2.append(data)
# Keeps corresponding numerical values. 
# SPs
c = dt.datetime.now()
current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')

i = 0
j = 0
k = 0
l = 0
SP_dict = {}
MT_dict = {}
CE_dict = {}

for Data in data_outputs_2:

    Data = h5.File(Data.filename, 'r')

    SP_d, MT_d, CE_d = utils.extract_masked_systems_by_seed(Data, survivors) 
    SP_dict.update(SP_d)
    MT_dict.update(MT_d)
    CE_dict.update(CE_d)

for Data in data_outputs_1:

    Data = h5.File(Data.filename, 'r')
    seeds = Data['SEED'][()]
    seed = seeds[0]  # Assuming we want to process the first seed in the file

    # Check if the seed exists in the survivors list
    if str(seed) not in survivors:
        print(f"Seed {seed} not in survivors, skipping...")
        i = i + 1
        continue
        fit_filename = pathToFiles + f"/stage_02__outputs_{seed[0]}.fits"

    fit_filename = pathToFiles + f"/stage_02_outputs_{seed}.fits" 
    hdu_pr = fits.PrimaryHDU()
    hdu_pr.writeto(fit_filename, overwrite=True)
    hdu = fits.open(fit_filename, mode='update')

    stellar_type_1 = Data['Stellar_Type(1)'][()]
    stellar_type_2 = Data['Stellar_Type(2)'][()]
    radius_1 = Data['Radius(1)'][()]
    radius_2 = Data['Radius(2)'][()]
    time = Data['Time'][()]
    semimajoraxis = (Data['SemiMajorAxis'][()]*const.R_sun).to(u.au).value
    roche_lobe_1 = Data['RocheLobe(1)'][()]
    roche_lobe_2 = Data['RocheLobe(2)'][()]
    mass_1 = Data['Mass(1)'][()]
    mass_2 = Data['Mass(2)'][()]
    mass_zams_1 = Data['Mass@ZAMS(1)'][()]
    mass_zams_2 = Data['Mass@ZAMS(2)'][()]
    mt_history = Data['MT_History'][()]
    eccentricity = Data['Eccentricity'][()]

    table_one = fits.BinTableHDU(Table(data=[seeds, time, stellar_type_1, stellar_type_2, semimajoraxis, radius_1, radius_2,
                                        roche_lobe_1, roche_lobe_2, mass_1, mass_2, mass_zams_1, mass_zams_2, mt_history, eccentricity],
                                    names=['Seed', 'Time', 'Stellar_Type_1', 'Stellar_Type_2', 'SemiMajorAxis',
                                            'Radius(1)', 'Radius(2)', 'RocheLobe(1)', 'RocheLobe(2)', 'Mass(1)', 'Mass(2)',
                                            'Mass@ZAMS(1)', 'Mass@ZAMS(2)', 'MT_History', 'Eccentricity'],
                                    units=[' ', 'Myr', ' ', ' ', 'AU', 'Rsun', 'Rsun', 'Rsun', 'Rsun', 'Msun', 'Msun', 'Msun', 'Msun', ' ', ' ']))
    hdu.append(table_one)
    hdu.flush()

    seedSP = [SP_dict[seed]['SEED'][()]]  # Convert to list for consistency
    statusSP = [SP_dict[seed]['Evolution_Status'][()]]
    stellarTypeZamsSP1 = [SP_dict[seed]['Stellar_Type@ZAMS(1)'][()]]
    stellarTypeZamsSP2 = [SP_dict[seed]['Stellar_Type@ZAMS(2)'][()]]
    stellarTypeSP1 = [SP_dict[seed]['Stellar_Type(1)'][()]]
    stellarTypeSP2 = [SP_dict[seed]['Stellar_Type(2)'][()]]
    massZamsSP1 = [SP_dict[seed]['Mass@ZAMS(1)'][()]]
    massZamsSP2 = [SP_dict[seed]['Mass@ZAMS(2)'][()]]
    semimajorAxisZamsSP = [SP_dict[seed]['SemiMajorAxis@ZAMS'][()]]  # in AU
    eccentricityZamsSP = [SP_dict[seed]['Eccentricity@ZAMS'][()]]
    mergersSP = [SP_dict[seed]['Merger'][()]]

    print(f'Lengths of SP arrays: {len(seedSP)}, {len(statusSP)}, {len(stellarTypeZamsSP1)}, {len(stellarTypeZamsSP2)}, '
          f'{len(stellarTypeSP1)}, {len(stellarTypeSP2)}, {len(massZamsSP1)}, {len(massZamsSP2)}, '
          f'{len(semimajorAxisZamsSP)}, {len(eccentricityZamsSP)}, {len(mergersSP)}')

    table_two = fits.BinTableHDU(Table(data=[seedSP, statusSP, stellarTypeZamsSP1, stellarTypeZamsSP2, stellarTypeSP1, stellarTypeSP2,
                                        massZamsSP1, massZamsSP2, semimajorAxisZamsSP, eccentricityZamsSP, mergersSP],
                                    names=['SEED', 'Evolution_Status', 'Stellar_Type@ZAMS(1)', 'Stellar_Type@ZAMS(2)',
                                            'Stellar_Type(1)', 'Stellar_Type(2)', 'Mass@ZAMS(1)', 'Mass@ZAMS(2)',
                                            'SemiMajorAxis@ZAMS', 'Eccentricity@ZAMS', 'Merger'],
                                    units=[' ', ' ', ' ', ' ', ' ', ' ', 'Msun', 'Msun', 'AU', '', '']))
    hdu.append(table_two)
    hdu.flush()

    if seed not in MT_dict:
        print(f"Seed {seed} not in MT_dict, third table not created.")
        j = j + 1
    
    else:

        seedMT = [MT_dict[seed]['SEED_MT'][()]]
        eventsMT = [MT_dict[seed]['MT_Event_Counter'][()]]
        stellarTypepreMT1 = [MT_dict[seed]['Stellar_Type(1)<MT'][()]]
        stellarTypepreMT2 = [MT_dict[seed]['Stellar_Type(2)<MT'][()]]
        stellarTypepstMT1 = [MT_dict[seed]['Stellar_Type(1)>MT'][()]]
        stellarTypepstMT2 = [MT_dict[seed]['Stellar_Type(2)>MT'][()]]
        masspreMT1 = [MT_dict[seed]['Mass(1)<MT'][()]]
        masspreMT2 = [MT_dict[seed]['Mass(2)<MT'][()]]
        masspstMT1 = [MT_dict[seed]['Mass(1)>MT'][()]]
        masspstMT2 = [MT_dict[seed]['Mass(2)>MT'][()]]
        semimajorAxispreMT = [MT_dict[seed]['SemiMajorAxis<MT'][()]]  # in AU
        semimajorAxispstMT = [MT_dict[seed]['SemiMajorAxis>MT'][()]]
        eccentricitypreMT = [MT_dict[seed]['Eccentricity<MT'][()]]
        eccentricitypstMT = [MT_dict[seed]['Eccentricity>MT'][()]]

        timepreMT = [MT_dict[seed]['Time<MT'][()]]
        timepstMT = [MT_dict[seed]['Time>MT'][()]]
        radiuspreMT1 = [MT_dict[seed]['Radius(1)<MT'][()]]
        radiuspreMT2 = [MT_dict[seed]['Radius(2)<MT'][()]]
        radiuspstMT1 = [MT_dict[seed]['Radius(1)>MT'][()]]
        radiuspstMT2 = [MT_dict[seed]['Radius(2)>MT'][()]]
        CEafterMT = [MT_dict[seed]['CE_History'][()]]

        table_three = fits.BinTableHDU(Table(data=[seedMT, eventsMT, stellarTypepreMT1, stellarTypepreMT2, stellarTypepstMT1, stellarTypepstMT2,
                                            masspreMT1, masspreMT2, masspstMT1, masspstMT2, semimajorAxispreMT, semimajorAxispstMT,
                                            eccentricitypreMT, eccentricitypstMT, timepreMT, timepstMT,
                                            radiuspreMT1, radiuspreMT2, radiuspstMT1, radiuspstMT2, CEafterMT],
                                        names=['SEED', 'MT_Event_Counter', 'Stellar_Type(1)<MT', 'Stellar_Type(2)<MT',
                                                'Stellar_Type(1)>MT', 'Stellar_Type(2)>MT', 'Mass(1)<MT', 'Mass(2)<MT',
                                                'Mass(1)>MT', 'Mass(2)>MT', 'SemiMajorAxis<MT', 'SemiMajorAxis>MT',
                                                'Eccentricity<MT', 'Eccentricity>MT', 'Time<MT', 'Time>MT',
                                                'Radius(1)<MT', 'Radius(2)<MT', 'Radius(1)>MT', 'Radius(2)>MT',
                                                'CE_History'],
                                        units=[' ', '', '', '', '', '', 'Msun', 'Msun', 'Msun', 'Msun', 'AU', 'AU',
                                                '', '', 'Myr', 'Myr', 'Rsun', 'Rsun', 'Rsun', 'Rsun', '']))
        hdu.append(table_three)
        hdu.flush()  

    if seed not in CE_dict:
        print(f"Seed {seed} not in CE_dict, fourth table not created.")
        k = k + 1
    else:

        seedCE = [CE_dict[seed]['SEED_CE'][()]]
        eventsCE = [CE_dict[seed]['CE_Event_Counter'][()]]
        stellarTypepreCE1 = [CE_dict[seed]['Stellar_Type(1)<CE'][()]]
        stellarTypepreCE2 = [CE_dict[seed]['Stellar_Type(2)<CE'][()]]
        stellarTypepstCE1 = [CE_dict[seed]['Stellar_Type(1)'][()]]
        stellarTypepstCE2 = [CE_dict[seed]['Stellar_Type(2)'][()]]
        masspreCE1 = [CE_dict[seed]['Mass(1)<CE'][()]]
        masspreCE2 = [CE_dict[seed]['Mass(2)<CE'][()]]
        masspstCE1 = [CE_dict[seed]['Mass(1)>CE'][()]]
        masspstCE2 = [CE_dict[seed]['Mass(2)>CE'][()]]
        semimajorAxispreCE = [CE_dict[seed]['SemiMajorAxis<CE'][()]]  # in au
        semimajorAxispstCE = [CE_dict[seed]['SemiMajorAxis>CE'][()]]
        timeCE = [CE_dict[seed]['Time_CE'][()]]
        circularizationTime = [CE_dict[seed]['Tau_Circ'][()]]
        radiuspreCE1 = [CE_dict[seed]['Radius(1)<CE'][()]]
        radiuspreCE2 = [CE_dict[seed]['Radius(2)<CE'][()]]
        radiuspstCE1 = [CE_dict[seed]['Radius(1)>CE'][()]]
        radiuspstCE2 = [CE_dict[seed]['Radius(2)>CE'][()]]

        table_four = fits.BinTableHDU(Table(data=[seedCE, eventsCE, stellarTypepreCE1, stellarTypepreCE2, stellarTypepstCE1, stellarTypepstCE2,
                                            masspreCE1, masspreCE2, masspstCE1, masspstCE2, semimajorAxispreCE, semimajorAxispstCE,
                                            timeCE, circularizationTime,
                                            radiuspreCE1, radiuspreCE2, radiuspstCE1, radiuspstCE2],
                                        names=['SEED', 'CE_Event_Counter', 'Stellar_Type(1)<CE', 'Stellar_Type(2)<CE',
                                                'Stellar_Type(1)>CE', 'Stellar_Type(2)>CE', 'Mass(1)<CE', 'Mass(2)<CE',
                                                'Mass(1)>CE', 'Mass(2)>CE', 'SemiMajorAxis<CE', 'SemiMajorAxis>CE',
                                                'Time_CE', 'Tau_Circ',
                                                'Radius(1)<CE', 'Radius(2)<CE', 'Radius(1)>CE', 'Radius(2)>CE'],
                                        units=[' ', '', '', '', '', '', 'Msun', 'Msun', 'Msun', 'Msun', 'AU', 'AU',
                                                'Myr', 'Myr', 'Rsun', 'Rsun', 'Rsun', 'Rsun']))
        hdu.append(table_four)
        hdu.flush()

    print(f"System of {l} (of {len(data_outputs_1)} sys.)")
    # Close the FITS file
    hdu.close()


    print(f"Data successfully written to {fit_filename}")

    l = l + 1
print(f"Total number of systems with BHs: {l}")
print(f"Total number of systems without MT data: {j}")
print(f"Total number of systems without CE data: {k}")
print(f"Total number of systems without CE data: {k}")


e = dt.datetime.now()
# Displays Time
current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
duration = e - s
print(" \n Finish time :", current_time)