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

matplotlib.use("Agg")


pathToData = '/data/a.saricaoglu/repo/COMPAS'

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
mode = ["Default_WD_Enabled/" ] #,"Limited_WD_Enabled/", "Default_WD_Disabled/", "Limited_WD_Disabled/"]    

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

    runs= [x[0] for x in os.walk(pathToData) if "Detailed_Output" in x[1]]
    print('runs:', runs)
    data_outputs = []
    for run in runs:
        out = [f for f in os.listdir(run + '/Detailed_Output') if ".h5" in f]
        for f in out:
            print("Reading file: ", run + '/Detailed_Output' + "/" + f)
            try:
                data = h5.File( run + '/Detailed_Output' + "/" + f)
                data_outputs.append(data)
            except:
                continue

    # Keeps corresponding numerical values. 
    # SPs
    if not os.path.exists(pathToData + '/Files/' + mod +  str(s.strftime("%m.%d"))): 
        os.makedirs(pathToData + '/Files/' + mod  +  str(s.strftime("%m.%d")))
    directoryf = pathToData + '/Files/' + mod  +  str(s.strftime("%m.%d"))

    i = 0

    c = dt.datetime.now()
    current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
    if not os.path.exists(pathToData+ "/Plots/" + mod +  str(c.strftime("%m.%d")+ "/" + current_time + "/") ): 
        os.makedirs(pathToData + "/Plots/"  + mod +  str(c.strftime("%m.%d")+ "/" + current_time + "/") ) 
    directoryp = pathToData + "/Plots/"  + mod +  str(c.strftime("%m.%d") + "/" + current_time + "/")  
    for Data in data_outputs:


        Data = h5.File(Data.filename, 'r')
        print("Batch (of " + str(len(data_outputs)) + " sys.) " + str(i) +  " start time :", current_time)

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

  
 