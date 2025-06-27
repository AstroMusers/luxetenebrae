#### PARENT SCRIPT: Stage_03_detailed_evolution.py
#----> this version plots individual properties.
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
from astropy.io import fits
from astropy.table import Table
from astropy import constants as const
from astropy import units as u
from luxetenebrae import calculations as calc   # functions from calculations.py
from luxetenebrae import utils as utils      # functions from utils.py

now = dt.datetime.now()
# Choose the modee to process
mode = 'WD_Enabled_Detailed'

# Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print("COMPAS_ROOT_DIR:", compasRootDir)
print(sys.path)

pathToData = os.path.join(compasRootDir, "luxetenebrae/runs", mode)
pathToFiles = os.path.join(compasRootDir, "luxetenebrae/files", mode)
pathToLogs = os.path.join(compasRootDir, "luxetenebrae/logs", mode)
pathToPlots = os.path.join(compasRootDir, "luxetenebrae/plots", mode)
# Create output directory if it doesn't exist
os.makedirs(pathToFiles, exist_ok=True)
os.makedirs(pathToLogs, exist_ok=True)
os.makedirs(pathToPlots, exist_ok=True)

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

out = [f for f in os.listdir(pathToFiles) if "stage_02_outputs" in f]
print('Output files found:', out)
detailed_output = []
SP_output = []
MT_output = []
CE_output = []
for f in out:
    print("Reading file: ", pathToFiles + "/" + f)

    data = fits.open(pathToFiles + "/" + f)
    print("Number of tables in the file:", len(data))
    detailed_output.append(data[1])
    SP_output.append(data[2])

    if len(data) > 3:
        MT_output.append(data[3])

    if len(data) > 4:
        CE_output.append(data[4])



# Keeps corresponding numerical values. 
# SPs
c = dt.datetime.now()
current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')

i = 0

merger_flags_manual = []
merger_flags_mthist = []
merger_flags = []
bhms_final_systems = []
bhbh_final_systems = []


for Data, Data_SP in zip(detailed_output, SP_output):

    Data = Data.data
    Data_SP = Data_SP.data
    # print("Batch (of " + str(len(data_outputs)) + " sys.) " + str(i) +  " start time :", current_time)
    seed = Data['Seed'][0]
    print('seed', seed)
    stellar_type_1 = Data['Stellar_Type_1']
    stellar_type_2 = Data['Stellar_Type_2']
    radius_1 = (Data['Radius(1)']*const.R_sun).to(u.au).value
    radius_2 = (Data['Radius(2)']*const.R_sun).to(u.au).value
    time = Data['Time']
    semimajoraxis = Data['SemiMajorAxis']
    roche_lobe_1 = (Data['RocheLobe(1)']*const.R_sun).to(u.au).value
    roche_lobe_2 = (Data['RocheLobe(2)']*const.R_sun).to(u.au).value
    mass_1 = Data['Mass(1)']
    mass_2 = Data['Mass(2)']
    mass_zams_1 = Data['Mass@ZAMS(1)']
    mass_zams_2 = Data['Mass@ZAMS(2)']
    mt_history = Data['MT_History']
    eccentricity = Data['Eccentricity']
    periapsiss = calc.periapsis(semimajoraxis, eccentricity)
    # rocheR_periapsis = calculate_periapsis()

    merger = Data_SP['Merger']

    mask = np.array([utils.is_ms_bh_pair(st1, st2) for st1, st2 in zip(stellar_type_1, stellar_type_2)])
    # print(f"Number of MS + BH states in this system: {np.sum(mask)}")

    merger_flag_manual = False
    for r1, r2, sa, e in zip(radius_1, radius_2, semimajoraxis, eccentricity):
        # print(check_merger_manual(r1, r2, sa))
        if utils.check_merger_manual(r1, r2, sa, e):
            print('check merger manual True')
            merger_flag_manual = True
    if merger_flag_manual:       
        merger_flags_manual.append(seed)

    merger_flag_mthist = False
    # print(check_merger_MT_hist(mt_history))
    if utils.check_merger_MT_hist(mt_history):
        print('check merger mt hist True')
        merger_flag_mthist = True
        print('r1 + r2 =', radius_1[-1] + radius_2[-1], 'semimajoraxis =', semimajoraxis[-1])
        merger_flags_mthist.append(seed)


    if utils.is_ms_bh_pair(stellar_type_1[-1], stellar_type_2[-1]):
        bhms_final_systems.append(seed)

    if utils.is_bh_bh_pair(stellar_type_1[-1], stellar_type_2[-1]):
        bhbh_final_systems.append(seed)

    merger_flags.append(merger)


    time_filtered = Data['Time'][:][mask]  # Assuming 'Time' is the time evolution data
    stype1_filtered = stellar_type_1[mask]
    stype2_filtered = stellar_type_2[mask]
    radius1_filtered = radius_1[mask]
    radius2_filtered = radius_2[mask]
    semimajoraxis_filtered = semimajoraxis[mask]
    bhms_lifetime = time_filtered[-1] - time_filtered[0]
    bhms_life = '{:.2f}'.format(bhms_lifetime) 
    semaj_average = np.mean(semimajoraxis_filtered)
    semaj_ave= '{:.2f}'.format(semaj_average) 

    fig, ax = plt.subplots(4,1, sharex=True, figsize=(4, 10))
    ax[0].scatter(time, stellar_type_1, s=2, label='Primary Star', color='blue')
    ax[0].scatter(time, stellar_type_2, s=2, label='Secondary Star', color='red')
    ax[0].set_yticks(np.arange(0,18,1), [r'MS (M<0.7 $M_\odot$)', r'MS (M>0.7 $M_\odot$)', 'HG', 'FGB', 'CHeB', 'EAGB', 'TPAGB', 'HeMS', 'HeHG', 'HeGB', 'HeWD', 'COWD', 'ONeWD','NS', 'BH', 'MR', 'CHE', 'None'])
    ax[0].fill_between(time_filtered, 0, np.max(stype1_filtered) , color='orange', alpha=0.1)
    ax[0].text(0.05, 0.95, fr'Lifetime of BH-MS: {bhms_life} Myr', transform=ax[0].transAxes,  verticalalignment='top')
    # ax[0].set_title('Stellar Type Evolution for MS + BH Systems')
    ax[0].legend()
    ax[0].grid(False) 
    ax[0].set_ylabel('Stellar Type')


    ax[1].plot(time, semimajoraxis, label='Semi-major Axis', color='yellow')
    ax[1].plot(time, radius_1, label=r'$R_p$', color='blue')
    ax[1].plot(time, radius_2, label=r'$R_s$', color='red')
    ax[1].plot(time, roche_lobe_1, label=r'$Roche R_p$', color='blue', linestyle='dashed')
    ax[1].plot(time, roche_lobe_2, label=r'$Roche R_s$', color='red', linestyle='dashed')   
    ax[1].plot(time, periapsiss, label='Periapsis', color='yellow', linestyle='dashed')

    ax[1].legend()
    ax[1].grid(False) 
    ax[1].fill_between(time_filtered, 0, max(np.max(semimajoraxis_filtered), np.max(radius_1), np.max(radius_2)), color='orange', alpha=0.1)
    ax[1].text(0.05, 0.15, rf'$\langle$SA$\rangle_\tau$: {semaj_ave} AU', transform=ax[1].transAxes,  verticalalignment='top')
    ax[1].set_ylabel('Semi-major Axis [AU]')

    ax[2].scatter(time, mt_history, s=2)
    ax[2].fill_between(time_filtered, 0, 7, color='orange', alpha=0.1)
    ax[2].set_yticks(np.arange(0,7,1), ['No MT', 'MT 1->2', 'MT 2->1', 'MTCE 1->2', 'MTCE 2->1', 'MTCE DoubleCore', 'Merger'])

    ax[2].text(0.05, 0.95, f'Merger (from flag): {merger}', transform=ax[2].transAxes,  verticalalignment='top')
    ax[2].text(0.05, 0.85, f'Merger (from mthist): {merger_flag_mthist}', transform=ax[2].transAxes,  verticalalignment='top')
    ax[2].text(0.05, 0.75, f'Merger (manual): {merger_flag_manual}', transform=ax[2].transAxes,  verticalalignment='top')
    # ax[2].legend()
    ax[2].grid(False) 
    ax[2].set_ylabel('MT History')

    ax[3].plot(time, eccentricity, label='Eccentricity', color='yellow')
    # ax[3].legend()
    ax[3].grid(False) 
    ax[3].set_ylabel('Eccentricity')        



    plt.xlabel('Time [Myr]')
    plt.legend()
    plt.grid(False) 
    plt.savefig(f'{pathToPlots}/stellar_type_evolution_{i}.png',bbox_inches='tight')
    plt.close()

    i = i + 1

print(f'Systems with merger (R1+R2>SA): {(merger_flags_manual)}')
print(f'Systems with merger (MT hist): {(merger_flags_mthist)}')
print(f'Systems with merger (flag): {(merger_flags)}')
print(f'Systems with BH-MS final state:{(bhms_final_systems)}')



print(f'Number of systems with merger (R1+R2>SA): {len(merger_flags_manual)}')
print(f'Number of systems with merger (MT hist): {len(merger_flags_mthist)}')
print(f'Number of systems with merger (flag): {np.sum(merger_flags)}')
print(f'Number of systems with BH-MS final state:{len(bhms_final_systems)}')
print(f'Number of systems with BH-BH final state:{len(bhbh_final_systems)}')
e = dt.datetime.now()
# Displays Time
current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
duration = e - s
print(" \n Finish time :", current_time)