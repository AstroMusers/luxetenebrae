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
from astropy.io import fits
from astropy.table import Table
from astropy import constants as const
from astropy import units as u
from luxetenebrae.utils import utils as utils      # functions from utils.py

now = dt.datetime.now()
# Choose the mode to process
mode = 'WD_Enabled_Detailed'

# Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
# print("COMPAS_ROOT_DIR:", compasRootDir)
# print(sys.path)

pathToData = os.path.join(compasRootDir, "luxetenebrae/runs", mode)
pathToLogs = os.path.join(compasRootDir, "luxetenebrae/logs", mode)
pathToFiles = os.path.join(compasRootDir, "luxetenebrae/files", mode)
# Create output directory if it doesn't exist
os.makedirs(pathToLogs, exist_ok=True)
os.makedirs(pathToFiles, exist_ok=True)
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

def purge_nonbh_systems():
    """    Purge non-black hole systems from the detailed output files."""
    # Displays Time

    s = dt.datetime.now()
    starttime = t.ctime(t.time())
    start = t.process_time()
    start_time = s.strftime("%d%m%y") + "_" + s.strftime('%H%M')
    print("Start time :", start_time)

    # Choose an output hdf5 file to work with
    runs= [x[0] for x in os.walk(pathToData) if "Detailed_Output" in x[1]]
    print('runs:', runs)
    data_outputs_1 = []
    for run in runs:  # Adjust the slice as needed
        out1 = [f for f in os.listdir(run + '/Detailed_Output') if ".h5" in f]
        for f in out1:
            print("Reading file: ", run + '/Detailed_Output' + "/" + f)
            try:
                data = h5.File(run + '/Detailed_Output' + "/" + f)
                data_outputs_1.append(data)
            except:
                continue

    # We save the seed number of the systems that survived the purge.
    survivors = open(f'{compasRootDir}/luxetenebrae/files/survivors.txt', 'a')

    i = 0
    j = 0
    k = 0
    l = 0

    # Loop through the data outputs and check for black holes and MS + BH pairs. If not found, delete the detailed output file. 
    print(f"Total number of systems to process: {len(data_outputs_1)}")
    print("Starting the purge...")

    for Data in data_outputs_1:

        k += 1
        Data = h5.File(Data.filename, 'r')
        # print(Data.keys())
        seed = Data['SEED'][()]

        print(f"System of {k} (of {len(data_outputs_1)} sys.)")

        stellar_type_1 = Data['Stellar_Type(1)'][:]
        stellar_type_2 = Data['Stellar_Type(2)'][:]
        BH1 = [x for x in stellar_type_1 if x == 14]
        BH2 = [x for x in stellar_type_2 if x == 14]

        print(f"System of {k} (of {len(data_outputs_1)} sys.) - Number of primary BH: {len(BH1)}, Number of secondary BH: {len(BH2)}")

        # Check if there are any BHs in the system
        if (len(BH1) + len(BH2)) == 0:
            print(f"No BHs in this system, deleting {Data.filename}")
            os.remove(f'{Data.filename}')
            i += 1
            continue

        # Check if there are any MS + BH pairs in the system
        mask = np.array([utils.is_ms_bh_pair(st1, st2) for st1, st2 in zip(stellar_type_1, stellar_type_2)])
        print(f"Number of MS + BH states in this system: {np.sum(mask)}")
        if np.sum(mask) == 0:
            print(f"No BH-MS states occur in this batch, deleting {Data.filename}")
            os.remove(f'{Data.filename}')
            j += 1
            continue
        
        # Save the seed of the system
        print(f"System {seed[0]} survived the purge.")
        survivors.write(f"{seed[0]}\n")
        l = l + 1
                
    survivors.close()

    print(f"Total number of systems: {k}")
    print(f"Total number of systems without BHs: {i}")
    print(f"Total number of systems without BH-MS pairs: {j}")           
    print(f"Total number of systems that survived the purge: {l}")
    
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    logging.shutdown()

    results = [f"Total number of systems: {k} Total number of systems without BHs: {i} Total number of systems without BH-MS pairs: {j} Total number of systems that survived the purge: {l} survivors file created at: {compasRootDir}/luxetenebrae/files/survivors.txt"
              ]

    return results
    
    