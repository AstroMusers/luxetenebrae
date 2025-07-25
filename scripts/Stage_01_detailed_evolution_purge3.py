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
from luxetenebrae import calculations as calc      # functions from calculations.py
from luxetenebrae.utils import is_ms_bh_pair

matplotlib.use("Agg")


pathToData = '/data/a.saricaoglu/repo/COMPAS/'

# Get the script name
script_name = os.path.basename(__file__)
# # Configure logging
# log_filename = f"{pathToData}/Files/{dt.datetime.now().strftime('%m.%d')}/{dt.datetime.now().strftime('%H%M')}/{script_name}_script.log"
# # os.makedirs(os.path.dirname(log_filename), exist_ok=True)
# logging.basicConfig(filename=log_filename, level=logging.INFO, 
#                     format='%(asctime)s - %(levelname)s - %(message)s')

# Redirect stdout and stderr to the log file
# class StreamToLogger:
#     def __init__(self, logger, log_level):
#         self.logger = logger
#         self.log_level = log_level
#         self.linebuf = ''

#     def write(self, buf):
#         for line in buf.rstrip().splitlines():
#             self.logger.log(self.log_level, line.rstrip())

#     def flush(self):
#         pass

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
    matplotlib.rcParams['figure.figsize'] = (15,10)
    matplotlib.rcParams['lines.markersize'] = 1
    matplotlib.rcParams['font.size'] = 14
    matplotlib.rcParams['legend.loc'] = "upper right"

    # runs= [x[0] for x in os.walk(pathToData) if "COMPAS_Output" in x[0]]
    runs = ['/data/a.saricaoglu/repo/COMPAS/COMPAS_Output_1']
    print('runs:', runs)
    data_outputs1 = []
    data_outputs2 = []
    for run in runs:
        out1 = [f for f in os.listdir(run + '/Detailed_Output') if ".h5" in f]
        out2 = [f for f in os.listdir(run) if ".h5" in f]
        for f in out1:
            print("Reading file: ", run + '/Detailed_Output' + "/" + f)
            try:
                data = h5.File( run + '/Detailed_Output' + "/" + f)
                data_outputs1.append(data)
            except:
                continue
                
        # for f in out2:
        #     print("Reading file: ", run  + "/" + f)
        #     data = h5.File( run + "/" + f)
        #     data_outputs2.append(data)
        #     # except:
        #     #     print("Error")
        #     #     continue


    i = 0

    c = dt.datetime.now()


    ps = open(f'{pathToData}Files/{mod}processed_systems.txt', 'r')    
    psr = ps.read().splitlines()

    # for Data in data_outputs2:
    #     Data = h5.File(Data.filename, 'r')
    #     print(Data.keys())
    #     SPs = Data['BSE System Paramteres'][()]      
    #     print(f'Size of the System:', len(SPs['SEED'])) 

    for Data in data_outputs1:


        Data = h5.File(Data.filename, 'r')
        # print(Data.keys())
        seed = Data['SEED'][()]

        # print(psr)
        if str(seed[0]) in psr:
            print(f"Seed {seed[0]} already processed, skipping...")
            continue
        print("Batch (of " + str(len(data_outputs1)) + " sys.) " + str(i) +  " start time :", current_time)

        stellar_type_1 = Data['Stellar_Type(1)'][:]
        stellar_type_2 = Data['Stellar_Type(2)'][:]

        semimajoraxis = (Data['SemiMajorAxis'][:]*const.R_sun).to(u.au).value
        BH1 = [x for x in stellar_type_1 if x == 14]
        BH2 = [x for x in stellar_type_2 if x == 14]
        print("BH1:", len(BH1), "BH2:", len(BH2))

        if (len(BH1) & len(BH2)) == 0:
            print(f"No BHs in this batch, deleting {Data.filename}")
            os.remove(f'{Data.filename}')
            continue

        mask = np.array([is_ms_bh_pair(st1, st2) for st1, st2 in zip(stellar_type_1, stellar_type_2)])
        print(f"Number of MS + BH states in this system: {np.sum(mask)}")
        if np.sum(mask) == 0:
            print(f"No BH-MS simultaneous emergence in this batch, deleting {Data.filename}")
            os.remove(f'{Data.filename}')
            continue



