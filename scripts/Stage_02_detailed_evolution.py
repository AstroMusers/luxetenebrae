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
import calculations as calc      # functions from calculations.py
from astropy.io import fits
from astropy.table import Table
from astropy import constants as const
from astropy import units as u
from utils import is_ms_bh_pair
# Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + 'postProcessing/PythonScripts')

matplotlib.use("Agg")


pathToData = '/data/a.saricaoglu/repo/COMPAS/'

# Get the script name
script_name = os.path.basename(__file__)
# Configure logging
log_filename = f"{pathToData}/Files/{dt.datetime.now().strftime('%m.%d')}/{dt.datetime.now().strftime('%H%M')}/{script_name}_script.log"
os.makedirs(os.path.dirname(log_filename), exist_ok=True)
logging.basicConfig(filename=log_filename, level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Redirect stdout and stderr to the log file
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

# sys.stdout = StreamToLogger(logging.getLogger('STDOUT'), logging.INFO)
# sys.stderr = StreamToLogger(logging.getLogger('STDERR'), logging.ERROR)
# Log the start of the script with the script name
logging.info(f'Script {script_name} started')
# Displays Time
s = dt.datetime.now()
starttime = t.ctime(t.time())
start = t.process_time()
start_time = s.strftime("%d%m%y") + "_" + s.strftime('%H%M')
print("Start time :", start_time)

# Choose the mode to process
mode = ["Default_WD_Enabled_Detailed/" ] #,"Limited_WD_Enabled/", "Default_WD_Disabled/", "Limited_WD_Disabled/"]    

# Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print(sys.path)

 

# Choose an output hdf5 file to work with
for mod in mode:
    # matplotlib.rcParams['figure.figsize'] = (15,10)
    # matplotlib.rcParams['lines.markersize'] = 1
    # matplotlib.rcParams['font.size'] = 14
    # matplotlib.rcParams['legend.loc'] = "upper right"

    ps = open(f'{pathToData}Files/{mod}survivors.txt', 'r')    
    survivors = ps.read().splitlines()

    
    runs= [x[0] for x in os.walk(pathToData + 'Runs/' + mod) if "Detailed_Output" in x[1]]
    print('runs:', runs)
    data_outputs = []
    data_outputs2 = []
    for run in runs:
        out = [f for f in os.listdir(run + '/Detailed_Output') if ".h5" in f]
        out2 = [f for f in os.listdir(run) if ".h5" in f]
        for f in out:
            print("Reading file: ", run + '/Detailed_Output' + "/" + f)
            try:
                data = h5.File( run + '/Detailed_Output' + "/" + f)
                data_outputs.append(data)
            except:
                continue
        for f in out2:
            print("Reading file: ", run  + "/" + f)
            data = h5.File( run + "/" + f, 'r')
            data_outputs2.append(data)
    # Keeps corresponding numerical values. 
    # SPs
    c = dt.datetime.now()
    current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
    if not os.path.exists(pathToData + '/Files/' + mod + '/' + str(s.strftime("%m.%d"))): 
        os.makedirs(pathToData + '/Files/' + mod + '/' +  str(s.strftime("%m.%d")))
    directoryf = pathToData + '/Files/' + mod  + '/' +   str(s.strftime("%m.%d"))


    


    # f = open(pathToData + '/Files/' + mod  + str(s.strftime("%m.%d")) +  "_general_outputs.txt", "a")

    # f.writelines(["\n","\n Run :", start_time])
    # print(f'Data outputs :, {data_outputs}')

    i = 0
    j = 0
    seed_mergerflag_dict = {}
    for Data in data_outputs2:
        Data = h5.File(Data.filename, 'r')
        SPs = Data['BSE_System_Parameters']   
        seedsSP = SPs['SEED'][()]
        mergersSP = SPs['Merger'][()]
        # If survivors mask is needed, apply it here
        seedsSP_mask = np.isin(seedsSP, survivors)  
        seeds_matched = seedsSP[seedsSP_mask]
        mergers_matched = mergersSP[seedsSP_mask]
        # Build the dictionary
        for seed, merger_flag in zip(seeds_matched, mergers_matched):
            seed_mergerflag_dict[str(seed)] = merger_flag  # Use str(seed) for consistent matching

        


    for Data in data_outputs:

        stellar_type_1 = Data['Stellar_Type(1)'][:]
        stellar_type_2 = Data['Stellar_Type(2)'][:]

        semimajoraxis = (Data['SemiMajorAxis'][:]*const.R_sun).to(u.au).value
        BH1 = [x for x in stellar_type_1 if x == 14]
        BH2 = [x for x in stellar_type_2 if x == 14]
        print("BH1:", len(BH1), "BH2:", len(BH2))

        if (len(BH1) & len(BH2)) == 0:
            print(f"No BHs in this batch, deleting {Data.filename}")
            os.remove(f'{Data.filename}')
            j += 1
            continue

        Data = h5.File(Data.filename, 'r')
        seed = Data['SEED'][()]
        seed_str = str(seed[0])  # Convert to string for matching
        merger_flag = seed_mergerflag_dict.get(seed_str, None)
        if merger_flag is not None:
            merger_flag_bool = bool(merger_flag)
            print(f"Seed: {seed_str}, Merger flag: {merger_flag}, As boolean: {merger_flag_bool}")

        merger = np.full(len(seed), merger_flag_bool)


        # print(psr)
        # if str(seed[0]) in psr:
        #     print(f"Seed {seed[0]} already processed, skipping...")
        #     continue
        ps = open(f'{pathToData}Files/{mod}survivors.txt', 'a')
        ps.write(f"{seed[0]}\n")
        ps.flush()
        ps.close()

        print("Batch (of " + str(len(data_outputs)) + " sys.) " + str(i) +  " start time :", current_time)

        
        # print(Data['Stellar_Type(1)'][()])
        # print(Data['Stellar_Type(2)'][()])
        # print(len(Data['Stellar_Type(1)'][()]))

        fit_filename = directoryf + f"/secundus_{seed[0]}.fits"
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

        table = fits.BinTableHDU(Table(data=[seed, time, stellar_type_1, stellar_type_2, semimajoraxis, radius_1, radius_2,
                                            roche_lobe_1, roche_lobe_2, mass_1, mass_2, mass_zams_1, mass_zams_2, mt_history, eccentricity, merger],
                                        names=['Seed', 'Time', 'Stellar_Type_1', 'Stellar_Type_2', 'SemiMajorAxis',
                                                'Radius(1)', 'Radius(2)', 'RocheLobe(1)', 'RocheLobe(2)', 'Mass(1)', 'Mass(2)',
                                                'Mass@ZAMS(1)', 'Mass@ZAMS(2)', 'MT_History', 'Eccentricity', 'Merger_Flag'], 
                                        units=[' ', 'Myr', ' ', ' ', 'AU', 'Rsun', 'Rsun', 'Rsun', 'Rsun', 'Msun', 'Msun', 'Msun', 'Msun', ' ', ' ', 'Boolean']))
        hdu.append(table)
        hdu.flush()
        # Close the FITS file
        hdu.close()


        print(f"Data successfully written to {fit_filename}")

    i = i + 1
print(f"Total number of systems with BHs: {i}")
print(f"Total number of systems without BHs: {j}")


e = dt.datetime.now()
# Displays Time
current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
duration = e - s
print(" \n Finish time :", current_time)