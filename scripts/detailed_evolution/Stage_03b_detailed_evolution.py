#### PARENT SCRIPT: Stage_03_detailed_evolution.py
#----> this version plots population properties.
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

semimajoraxis_average = []
lifetime_average = []
eccentricity_average = []

j = 0
k = 0

primary_bh = 0
secondary_bh = 0
# Choose an output hdf5 file to work with

ps = open(f'{pathToData}Files/{mod}survivors.txt', 'r')    
survivors = ps.read().splitlines()


out = [f for f in os.listdir(pathToData + '/Files/' + mod + day) if ".fits" in f]
data_outputs = []
for f in out:
    # print("Reading file: ", pathToData + '/Files/' + mod + day + "/" + f)

    data = fits.open(pathToData + '/Files/' + mod + day + "/" + f)
    data_outputs.append(data[1])

runs= [x[0] for x in os.walk(pathToData + 'Runs/' + mod) if "Detailed_Output" in x[1]]
print('runs:', runs)
data_outputs2 = []
for run in runs:
    out2 = [f for f in os.listdir(run) if ".h5" in f]
    for f in out2:
        print("Reading file: ", run  + "/" + f)
        data = h5.File( run + "/" + f, 'r')
        data_outputs2.append(data)

# Keeps corresponding numerical values. 
# SPs
c = dt.datetime.now()
current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')





# f = open(pathToData + '/Files/' + mod  + str(s.strftime("%m.%d")) +  "_general_outputs.txt", "a")

# f.writelines(["\n","\n Run :", start_time])
# print(f'Data outputs :, {data_outputs}')
primary_mass = []
secondary_mass = []
survivor_primary_mass = []
survivor_secondary_mass = []
for Data in data_outputs2:
    Data = h5.File(Data.filename, 'r')
    SPs = Data['BSE_System_Parameters']   
    massZamsSP1 = SPs['Mass@ZAMS(1)'][()] 
    massZamsSP2 = SPs['Mass@ZAMS(2)'][()]
    primary_mass = np.concatenate((primary_mass, massZamsSP1))
    secondary_mass = np.concatenate((secondary_mass, massZamsSP2))
    seedsSP = SPs['SEED'][()]
    seedsSP_mask = np.isin(seedsSP, survivors)  
    seeds_matched = seedsSP[seedsSP_mask]  
    survivor_primary_mass = np.concatenate((massZamsSP1[seedsSP_mask], survivor_primary_mass))
    survivor_secondary_mass = np.concatenate((massZamsSP2[seedsSP_mask], survivor_secondary_mass))     
    # Build the dictionary

for Dat in data_outputs:
    
    Data = Dat.data
    print("Batch (of " + str(len(data_outputs)) + " sys.) " + str(i) +  " start time :", current_time)
    semimajoraxis = Data['SemiMajorAxis']
    # print(f'Length of semimajor axis array: {len(semimajoraxis)}')
    # print(f'Minimum semimajor axis: {np.min(semimajoraxis)}')

    if np.any(semimajoraxis < 0):
        print(f"Negative semimajor axis found in batch {i}, skipping this batch")
        k += 1
        continue
    seed = Data['Seed'][0]
    print('seed', seed)
    stellar_type_1 = Data['Stellar_Type_1']
    stellar_type_2 = Data['Stellar_Type_2']
    radius_1 = (Data['Radius(1)']*const.R_sun).to(u.au).value
    radius_2 = (Data['Radius(2)']*const.R_sun).to(u.au).value
    time = Data['Time']
    roche_lobe_1 = (Data['RocheLobe(1)']*const.R_sun).to(u.au).value
    roche_lobe_2 = (Data['RocheLobe(2)']*const.R_sun).to(u.au).value
    mass_1 = Data['Mass(1)']
    mass_2 = Data['Mass(2)']
    mass_zams_1 = Data['Mass@ZAMS(1)']
    mass_zams_2 = Data['Mass@ZAMS(2)']
    mt_history = Data['MT_History']
    merger = Data['Merger_Flag'][0]
    eccentricity = Data['Eccentricity']
    periapsiss = periapsis(semimajoraxis, eccentricity)
    # rocheR_periapsis = calculate_periapsis()

    bhms_mask = np.array([is_ms_bh_pair(st1, st2) for st1, st2 in zip(stellar_type_1, stellar_type_2)])

    print(f"Number of MS + BH states in this system: {np.sum(bhms_mask)}")
    if np.sum(bhms_mask) == 0:
        print(f"No BH-MS simultaneous emergence in this batch, skipping")
        j += 1
        continue

    if is_ms_bh_pair(stellar_type_1[-1], stellar_type_2[-1]):
        bhms_final_systems.append(seed)

    if is_bh_bh_pair(stellar_type_1[-1], stellar_type_2[-1]):
        bhbh_final_systems.append(seed)            

    time_filtered = time[bhms_mask]   # Assuming 'Time' is the time evolution data
    stype1_filtered = stellar_type_1[bhms_mask] 
    stype2_filtered = stellar_type_2[bhms_mask] 
    radius1_filtered = radius_1[bhms_mask] 
    radius2_filtered = radius_2[bhms_mask] 
    semimajoraxis_filtered = semimajoraxis[bhms_mask] 
    eccentricity_filtered = eccentricity[bhms_mask] 

    bhms_lifetime = time_filtered[-1] - time_filtered[0]
    lifetime_average.append(bhms_lifetime) 

    semaj_average = np.mean(semimajoraxis_filtered)
    semimajoraxis_average.append(semaj_average)

    ecc_average = np.mean(eccentricity_filtered)
    eccentricity_average.append(ecc_average)

    primary_mask = np.array([which_bh(st1, st2) for st1, st2 in zip(stype1_filtered, stype2_filtered)])
    # primary_mass_filtered = mass_zams_1[primary_mask == 1]
    # secondary_mass_filtered = mass_zams_2[primary_mask == 0]
    if np.sum(primary_mask) < len(primary_mask):
        secondary_bh += 1
    else:
        primary_bh += 1

    i = i + 1

    fig, ax = plt.subplots(3,1, figsize=(3.5, 7))
    ax[0].hist(semimajoraxis_average, bins=50, color='blue', alpha=0.7, label='Semimajor Axis')
    ax[0].set_xlabel('Semimajor Axis (AU)')
    ax[0].set_ylabel('Number of Systems')
    ax[0].set_yscale('log')
    ax[1].hist(lifetime_average, bins=50, color='green', alpha=0.7, label='Lifetime')
    ax[1].set_xlabel('Lifetime (Myr)')
    ax[1].set_ylabel('Number of Systems')
    ax[1].set_yscale('log')
    ax[2].hist(eccentricity_average, bins=50, color='red', alpha=0.7, label='Eccentricity')
    ax[2].set_xlabel('Eccentricity')
    ax[2].set_ylabel('Number of Systems')
    ax[2].set_yscale('log')
    #plt.suptitle('BH-MS Systems Average Histogram', fontsize=16)
    plt.tight_layout()
    plt.savefig(directoryp + 'bhms_properties_histogram.png')
    plt.close(fig)


    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    # counts, bin_edges1 = np.histogram(primary_mass, bins=50)
    # counts, bin_edges2 = np.histogram(secondary_mass, bins=50)
    # counts1, bin_edges = np.histogram(survivor_primary_mass, bins=50)
    # counts2, bin_edges = np.histogram(survivor_secondary_mass, bins=50)
    # survivor_primary_mass_density = density(counts1, bins=bin_edges1)
    # survivor_secondary_mass_density = density(counts2, bins=bin_edges2)
    ax.hist(primary_mass, bins=50, color='blue', alpha=0.7, label='Primary Mass', histtype='step')
    ax.hist(survivor_primary_mass, bins=50, color='blue', alpha=0.7, label='BHMS Primary Mass', histtype='step', linestyle='--')
    ax.hist(secondary_mass, bins=50, color='orange', alpha=0.7, label='Secondary Mass', histtype='step')
    ax.hist(survivor_secondary_mass, bins=50, color='orange', alpha=0.7, label='BHMS Secondary Mass', histtype='step', linestyle='--')
    ax.set_xlabel('Mass (M☉)')
    ax.set_ylabel('Number of Systems')
    ax.set_yscale('log')
    ax.legend()
    #ax.set_title('Initial Mass Distribution of BH-MS Systems')
    plt.tight_layout()
    plt.savefig(directoryp + 'bhms_mass_distribution_initial.png')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    ax.hist(primary_mass, bins=50, color='blue', alpha=0.7, label='Primary Mass', density=True)
    ax.hist(secondary_mass, bins=50, color='orange', alpha=0.7, label='Secondary Mass', density=True)
    kroupa_fit = kroupa_imf_normalized(5, 150)
    ax.plot(np.linspace(5, 150, 1000), kroupa_fit, color='black', linestyle='--', label='Kroupa IMF')
    ax.set_xlabel('Mass (M☉)')
    ax.set_ylabel('Probability Density')
    ax.set_yscale('log')
    ax.legend()
    plt.tight_layout()
    plt.savefig(directoryp + 'bhms_mass_distribution_density.png')
    plt.close()

    print(f'Systems with merger (R1+R2>SA): {(merger_flags_manual)}')
    print(f'Systems with merger (MT hist): {(merger_flags_mthist)}')
    print(f'Systems with merger (flag): {(merger_flags)}')
    print(f'Systems with BH-MS final state:{(bhms_final_systems)}')



    print(f'Number of systems with merger (R1+R2>SA): {len(merger_flags_manual)}')
    print(f'Number of systems with merger (MT hist): {len(merger_flags_mthist)}')
    print(f'Number of systems with merger (flag): {np.sum(merger_flags)}')
    print(f'Number of systems with BH-MS final state:{len(bhms_final_systems)}')
    print(f'Number of systems with BH-BH final state:{len(bhbh_final_systems)}')

    print(f'Number of systems with primary BH: {primary_bh}')
    print(f'Number of systems with secondary BH: {secondary_bh}')
    print(f'Total number of systems with BHs: {i}')
    print(f'Total number of systems with BH-MS pairs: {len(lifetime_average)}')

    print(f'Skipped systems due to absence of BH-MS pairs: {j}')
    print(f'Skipped systems due to negative semimajor axis: {k}')

    print(f'Minimum primary mass: {np.min(survivor_primary_mass)} M☉')
    print(f'Maximum primary mass: {np.max(survivor_primary_mass)} M☉')
# Displays Time
e = dt.datetime.now()
current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
duration = e - s
print(" \n Finish time :", current_time)