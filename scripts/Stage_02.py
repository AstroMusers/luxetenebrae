#### PARENT SCRIPT: tertius_masking_020525
#----> this version, instead of the systems just before the first and last mass transfer (MT) event, marks MT events in which StellarType1 or StellarType2
#      and gets the system parameters before and after each event. CAUTION: since each transition event for the system is stored as another (with the same seed number)
#      seed, arrays will be larger than actual system number.
import os, sys
import numpy as np               # for handling arrays
import numpy.ma as ma
import pandas as pd
import h5py as h5                # for reading the COMPAS data
import time as t                      # for finding computation time
from datetime import datetime
import matplotlib.pyplot as plt  # for plotting
import matplotlib
from luxetenebrae import calculations as calc      # functions from calculations.py
from luxetenebrae import utils as utils
from astropy.io import fits
from astropy import units as u
from astropy import constants as const
from astropy.table import Table
import logging
import datetime as dt

now = dt.datetime.now()
# Choose the mode to process
mode = 'WD_Enabled'

# Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print("COMPAS_ROOT_DIR:", compasRootDir)
print(sys.path)


pathToFiles = os.path.join(compasRootDir, "luxetenebrae/files", mode)
pathToLogs = os.path.join(compasRootDir, "luxetenebrae/logs", mode)
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

fit_filename = pathToFiles + "/stage_02_output.fits"
hdu_pr = fits.PrimaryHDU()
hdu_pr.writeto(fit_filename, overwrite=True)
hdu = fits.open(fit_filename, mode='update')

# Boolen arrays for MASKing.
MASKSPBH1 = np.array([], dtype='bool')
MASKSPBH2 = np.array([], dtype='bool')

MASKPREMT = np.array([], dtype='bool')
MASKPSTMT = np.array([], dtype='bool')

MASKSPunb = np.array([], dtype='bool')
MASKSPdco = np.array([], dtype='bool')
MASKSPmrgr = np.array([], dtype='bool')

MASKPREMTBH1 = np.array([], dtype='bool')
MASKPSTMTBH1 = np.array([], dtype='bool')

MASKPREMTBH2 = np.array([], dtype='bool')
MASKPSTMTBH2 = np.array([], dtype='bool')

MASKPREMTinSP = np.array([], dtype='bool')
MASKPSTMTinSP = np.array([], dtype='bool')

MASKPREMTBHNS1 = np.array([], dtype='bool')
MASKPSTMTBHNS1 = np.array([], dtype='bool')

MASKPREMTBHNS2 = np.array([], dtype='bool')
MASKPSTMTBHNS2 = np.array([], dtype='bool')

MASKPREMTnonBH = np.array([], dtype='bool')
MASKPSTMTnonBH = np.array([], dtype='bool')

MASKPREMTSEMAJ = np.array([], dtype='bool')
MASKPSTMTSEMAJ = np.array([], dtype='bool')

MASKPREMTORBPER = np.array([], dtype='bool')
MASKPSTMTORBPER = np.array([], dtype='bool')

MASKPREMTnegCE = np.array([], dtype='bool')
MASKPSTMTnegCE = np.array([], dtype='bool')

MASKPREMTposCE = np.array([], dtype='bool')
MASKPSTMTposCE = np.array([], dtype='bool')

MASKFIRSTMT = np.array([], dtype='bool')
MASKLASTMT = np.array([], dtype='bool')

MASKPREMTSEARCHABILITY = np.array([], dtype='bool')
MASKPSTMTSEARCHABILITY = np.array([], dtype='bool')

fits_data = pathToFiles + '/stage_01_output.fits'
with fits.open(fits_data) as hdul:
    hdul.info()
    SP = hdul[1].data
    MT = hdul[2].data
    CE = hdul[3].data


# for file in files:
#     print("file :", file)
#     for var in filenames:
#         print("var :", var)
#         print("if check :" , arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]])
#         if var in file and len(arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]]) == 0 :
#             index = np.where(np.asarray(filenames) == var)[0][0]
#             arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]] = 1

#     print(pathToData  + file)
#     data = np.loadtxt(pathToData  + file)
#     print("index: ", index)
#     print("File name: ", file, " Variable name: ", filenames[index])
#     if len(data) != 0:
#         print("data reached.")
#     print("len: ", len(data))
#     print("df col name:", filenames[index])   

#     if index == 0:
#         lenSP = len(data)
#     if index == 1:
#         lenMT = len(data)
#     if index == 2:
#         lenCE = len(data_CE)

#     if len(data) == lenSP:

#         SPdf[filenames[index]] = data
#         print(SPdf[filenames[index]] .isnull().sum(), SPdf[filenames[index]] .max())
#         print("SUCCESS (SP)")
#         i = i+1

#     if len(data) == lenMT:

#         MTdf[filenames[index]] = data
#         print(MTdf[filenames[index]].isnull().sum(), MTdf[filenames[index]].max())
#         print("SUCCESS (MT)")
#         i = i+1

#     if len(data) == lenCE:

#         CEdf[filenames[index]] = data
#         print(CEdf[filenames[index]].isnull().sum(), CEdf[filenames[index]].max())
#         print("SUCCESS (CE)")
#         i = i+1


# Chooses first MT events.
MASKFIRSTMT = (MT['EVENTSMT'] == 1)

# Chooses last MT events.
MASKLASTMT = utils.find_last_mt(MT['MTs'])

# Chooses MT events where stellar type changes.
MASKTYPECHANGE1, MASKTYPECHANGE2 = utils.type_change_MT(MT['STELLARTYPEPREMT1'],MT['STELLARTYPEPSTMT1'],MT['STELLARTYPEPREMT2'], MT['STELLARTYPEPSTMT2'], start_time)


ORBITALPERIODSP = calc.orbital_period(SP['MASSZAMSSP1'], SP['MASSZAMSSP2'], SP['SEMIMAJORAXISZAMSSP'])
ORBITALPERIODPREMT = calc.orbital_period(MT['MASSPREMT1'], MT['MASSPREMT2'], MT['SEMIMAJORAXISPREMT'])
ORBITALPERIODPSTMT = calc.orbital_period(MT['MASSPSTMT1'], MT['MASSPSTMT2'], MT['SEMIMAJORAXISPSTMT'])
print(np.shape(ORBITALPERIODPREMT))

COSIPREMT = calc.orbital_inclination(MT['RADIUSPREMT1'], MT['RADIUSPREMT2'], MT['SEMIMAJORAXISPREMT'])
COSIPSTMT = calc.orbital_inclination(MT['RADIUSPSTMT1'], MT['RADIUSPSTMT2'], MT['SEMIMAJORAXISPSTMT'])
print(np.shape(COSIPREMT))
MASKPREMTSEARCHABILITY = calc.searchability(COSIPREMT)
MASKPSTMTSEARCHABILITY = calc.searchability(COSIPSTMT)

print(MASKPREMTSEARCHABILITY)
print(np.shape(MASKPREMTSEARCHABILITY))
# i = np.random.uniform(0,180,len(MTs))



## may add eccentricity later on
CEAFTERMT = MT['CEAFTERMT']

MASKSPBH1 = ((SP['STELLARTYPESP1']  == 14) & (SP['STELLARTYPESP2']  != 14)) 
MASKSPBH2 = ((SP['STELLARTYPESP1']  != 14) & (SP['STELLARTYPESP2']  == 14))
MASKPREMTBH1 = ((MT['STELLARTYPEPREMT1']  == 14) & (MT['STELLARTYPEPREMT2']  != 14)) 
MASKPREMTBH2 = ((MT['STELLARTYPEPREMT1']  != 14) & (MT['STELLARTYPEPREMT2']  == 14))       
MASKPSTMTBH1 = ((MT['STELLARTYPEPSTMT1']  == 14) & (MT['STELLARTYPEPSTMT2']  != 14)) 
MASKPSTMTBH2 = ((MT['STELLARTYPEPSTMT1']  != 14) & (MT['STELLARTYPEPSTMT2']  == 14))

MASKSPunb = (SP['EVOLUTIONSTATSP'] == 14)
MASKSPdco = (SP['EVOLUTIONSTATSP'] == 11)
MASKSPmrgr = (SP['EVOLUTIONSTATSP'] == 9) 

MASKMTinSP =  np.in1d(SP['SPs'], MT['MTs'])

MASKPREMTnonBH = ((MT['STELLARTYPEPREMT1'] != 14) & (MT['STELLARTYPEPREMT2'] != 14)) #none is bh
MASKPREMTBHNS1 = ((MT['STELLARTYPEPREMT1'] == 14) & (MT['STELLARTYPEPREMT2'] == 13)) #primary is bh secondary is ns
MASKPREMTBHNS2 = ((MT['STELLARTYPEPREMT1'] == 13) & (MT['STELLARTYPEPREMT2'] == 14)) #primary is ns secondary is bh
MASKPREMTBHMS1 = ((MT['STELLARTYPEPREMT1'] == 14) & ((MT['STELLARTYPEPREMT2'] == 0)|(MT['STELLARTYPEPREMT2'] == 1))) #primary is bh secondary is MS
MASKPREMTBHMS2 = (((MT['STELLARTYPEPREMT1'] == 0)|(MT['STELLARTYPEPREMT1'] == 1)) & (MT['STELLARTYPEPREMT2'] == 14)) #primary is MS secondary is bh

MASKPSTMTnonBH = ((MT['STELLARTYPEPSTMT1'] != 14) & (MT['STELLARTYPEPSTMT2'] != 14)) #none is bh
MASKPSTMTBHNS1 = ((MT['STELLARTYPEPSTMT1'] == 14) & (MT['STELLARTYPEPSTMT2'] == 13)) #primary is bh secondary is ns
MASKPSTMTBHNS2 = ((MT['STELLARTYPEPSTMT1'] == 13) & (MT['STELLARTYPEPSTMT2'] == 14)) #primary is ns secondary is bh
MASKPSTMTBHMS1 = ((MT['STELLARTYPEPSTMT1'] == 14) & ((MT['STELLARTYPEPSTMT2'] == 0) | (MT['STELLARTYPEPSTMT2'] == 1))) #primary is bh secondary is MS
MASKPSTMTBHMS2 = (((MT['STELLARTYPEPSTMT1'] == 0) | (MT['STELLARTYPEPSTMT1'] == 1)) & (MT['STELLARTYPEPSTMT2'] == 14)) #primary is MS secondary is bh


MASKPREMTSEMAJ = (MT['SEMIMAJORAXISPREMT'] <= 3.0)
MASKPSTMTSEMAJ = (MT['SEMIMAJORAXISPSTMT'] <= 3.0)

MASKPREMTORBPER = (ORBITALPERIODPREMT <= 70)
MASKPSTMTORBPER = (ORBITALPERIODPSTMT <= 70)

MASKPREMTSEMAJ_masuda = (MT['SEMIMAJORAXISPREMT'] <= 0.4)
MASKPSTMTSEMAJ_masuda = (MT['SEMIMAJORAXISPSTMT'] <= 0.4)

MASKPREMTORBPER_masuda = (ORBITALPERIODPREMT <= 30)
MASKPSTMTORBPER_masuda = (ORBITALPERIODPSTMT <= 30)

MASKPREMTnegCE =  np.in1d(MT['MTs']*MASKFIRSTMT, CE['CEs'], invert=True)
MASKPREMTposCE =  np.in1d(MT['MTs']*MASKFIRSTMT, CE['CEs'])
MASKPSTMTnegCE =  np.in1d(MT['MTs']*MASKLASTMT, CE['CEs'], invert=True)
MASKPSTMTposCE =  np.in1d(MT['MTs']*MASKLASTMT, CE['CEs'])
# print('maskpremtnegce, posce, pstnegce, pstposce ', len(MASKPREMTnegCE),len(MASKPREMTposCE),len(MASKPSTMTnegCE),len(MASKPSTMTposCE))

# arrays_to_save = [MASKPREMTBH1,MASKPREMTBH2,MASKPREMTBHNS1,MASKPREMTBHNS2,
#                 MASKPREMTnonBH, MASKPREMTSEMAJ, MASKPREMTnegCE, MASKPREMTposCE,MASKFIRSTMT,
#                 MASKPSTMTBH1,MASKPSTMTBH2,MASKPSTMTBHNS1,MASKPSTMTBHNS2, MASKPSTMTnonBH, MASKPSTMTSEMAJ,
#                 MASKPSTMTnegCE, MASKPSTMTposCE, MASKLASTMT, MASKTYPECHANGE1, MASKTYPECHANGE2,
#                 ORBITALPERIODPREMT,ORBITALPERIODPSTMT, MASKPREMTORBPER, COSIPREMT,MASKPSTMTORBPER, COSIPSTMT,
#                 MASKPREMTSEARCHABILITY, MASKPSTMTSEARCHABILITY, CEAFTERMT]

# filenames = [   "MASKPREMTBH1", "MASKPREMTBH2", "MASKPREMTBHNS1", "MASKPREMTBHNS2",
#                 "MASKPREMTnonBH", "MASKPREMTSEMAJ", "MASKPREMTnegCE", "MASKPREMTposCE", "MASKFIRSTMT",
#                 "MASKPSTMTBH1", "MASKPSTMTBH2", "MASKPSTMTBHNS1", "MASKPSTMTBHNS2", "MASKPSTMTnonBH", "MASKPSTMTSEMAJ", 
#                 "MASKPSTMTnegCE", "MASKPSTMTposCE", "MASKLASTMT", "MASKTYPECHANGE1", "MASKTYPECHANGE2",
#                 "ORBITALPERIODPREMT", "ORBITALPERIODPSTMT","MASKPREMTORBPER", "COSIPREMT", "MASKPSTMTORBPER", "COSIPSTMT",
#                 "MASKPREMTSEARCHABILITY", "MASKPSTMTSEARCHABILITY","CEAFTERMT"]

# for e in arrays_to_save:
#     print(len(e))


print(f'Run : {start_time}')

print(f"Number of systems: {len(SP['MASSZAMSSP1'])}") 
print(f"Number of systems forming MT (via SP STATUS): {np.sum(MASKMTinSP)}") 
print(f"Number of systems forming MT (via MT mass): {np.sum(MASKFIRSTMT)}")  
print(f"Number of primary mass BHs (SP):  {np.sum(MASKSPBH1)}") 
print(f"Number of primary mass BHs (MT): {np.sum(MASKPREMTBH1)}")   
print(f"Number of secondary mass BHs (SP): {np.sum(MASKSPBH2)}") 
print(f"Number of secondary mass BHs (MT): {np.sum(MASKPREMTBH2)}")      
print(f"Number of unbound binaries (SP): {np.sum(MASKSPunb)}")  
print(f"Number of double compact objects (SP): {np.sum(MASKSPdco)}")   
print(f"Number of mergers (SP): {np.sum(MASKSPmrgr)}")
print(f"Number of BH-NS binaries (pre MT): {np.sum(MASKPREMTBHNS1)}") 
print(f"Number of NS-BH binaries (pre MT): {np.sum(MASKPREMTBHNS2)}") 

print(f"Number of systems without any BHs (pre MT):  {np.sum(MASKPREMTnonBH)}") 
print(f"Number of systems with semimajor axis < 3.0 AU (pre MT):  {np.sum(MASKPREMTSEMAJ)}") 
print(f"Number of systems undergoing CE (pre MT):  {np.sum(MASKPREMTposCE)}") 
print(f"Number of systems not undergoing CE (pre MT):  {np.sum(MASKPREMTnegCE)}") 
print(f"Number of systems with semimajor axis < 3.0 AU AND undergoing CE (pre MT):  {np.sum(MASKPREMTSEMAJ*MASKPREMTposCE)}") 
print(f"Number of systems with semimajor axis < 3.0 AU AND NOT undergoing CE (pre MT):  {np.sum(MASKPREMTSEMAJ*MASKPREMTnegCE)}") 
print(f"Number of systems with orbital period < 70 days AND undergoing CE (pre MT):  {np.sum(MASKPREMTORBPER*MASKPREMTposCE)}") 
print(f"Number of systems with orbital period < 70 days AND NOT undergoing CE (pre MT):  {np.sum(MASKPREMTORBPER*MASKPREMTnegCE)}")     
print(f"Number of systems with orbital inclination (COSI) <  np.cos(np.deg2rad(89.9)) AND undergoing CE (pre MT):  {np.sum(MASKPREMTSEARCHABILITY*MASKPREMTposCE)}") 
print(f"Number of systems with orbital inclination (COSI) <  np.cos(np.deg2rad(89.9)) AND NOT undergoing CE (pre MT):  {np.sum(MASKPREMTSEARCHABILITY*MASKPREMTnegCE)}") 
print(f"Number of BH-NS binaries (pst MT):  {np.sum(MASKPSTMTBHNS1)}") 
print(f"Number of NS-BH binaries (pst MT):  {np.sum(MASKPSTMTBHNS2)}") 

print(f"Number of systems without any BHs (pst MT):  {np.sum(MASKPSTMTnonBH)}") 
print(f"Number of systems with semimajor axis < 3.0 AU (pst MT):  {np.sum(MASKPSTMTSEMAJ)}") 
print(f"Number of systems undergoing CE (pst MT):  {np.sum(MASKPSTMTposCE)}") 
print(f"Number of systems not undergoing CE (pst MT):  {np.sum(MASKPSTMTnegCE)}") 
print(f"Number of systems with semimajor axis < 3.0 AU AND undergoing CE (pst MT):  {np.sum(MASKPSTMTSEMAJ*MASKPSTMTposCE)}") 
print(f"Number of systems with semimajor axis < 3.0 AU AND NOT undergoing CE (pst MT):  {np.sum(MASKPSTMTSEMAJ*MASKPSTMTnegCE)}") 
print(f"Number of systems with orbital period < 70 days AND undergoing CE (pst MT):  {np.sum(MASKPSTMTORBPER*MASKPSTMTposCE)}") 
print(f"Number of systems with orbital period < 70 days AND NOT undergoing CE (pst MT):  {np.sum(MASKPSTMTORBPER*MASKPSTMTnegCE)}") 
print(f"Number of systems with orbital inclination (COSI) <  np.cos(np.deg2rad(89.9)) AND undergoing CE (pst MT):  {np.sum(MASKPSTMTSEARCHABILITY*MASKPSTMTposCE)}") 
print(f"Number of systems with orbital inclination (COSI) <  np.cos(np.deg2rad(89.9)) AND NOT undergoing CE (pst MT):  {np.sum(MASKPSTMTSEARCHABILITY*MASKPSTMTnegCE)}") 
print(f"Total MT events:  {np.sum(MT['EVENTSMT'])}") 
print(f"Number of primary HeWD (pre MT) :  {np.sum([MT['STELLARTYPEPREMT1'] == 10 ])}") 
print(f"Number of primary COWD (pre MT) :  {np.sum([MT['STELLARTYPEPREMT1'] == 11 ])}") 
print(f"Number of primary ONeWD (pre MT) :  {np.sum([MT['STELLARTYPEPREMT1'] == 12 ])}") 
print(f"Number of primary HeWD (pst MT) :  {np.sum([MT['STELLARTYPEPSTMT1'] == 10 ])}") 
print(f"Number of primary COWD (pst MT) :  {np.sum([MT['STELLARTYPEPSTMT1'] == 11 ])}") 
print(f"Number of primary ONeWD (pst MT) :  {np.sum([MT['STELLARTYPEPSTMT1'] == 12 ])}") 
print(f"Number of secondary HeWD (pre MT) :  {np.sum([MT['STELLARTYPEPREMT2'] == 10 ])}") 
print(f"Number of secondary COWD (pre MT) :  {np.sum([MT['STELLARTYPEPREMT2'] == 11 ])}") 
print(f"Number of secondary ONeWD (pre MT) :  {np.sum([MT['STELLARTYPEPREMT2'] == 12 ])}") 
print(f"Number of secondary HeWD (pst MT) :  {np.sum([MT['STELLARTYPEPSTMT2'] == 10 ])}") 
print(f"Number of secondary COWD (pst MT) :  {np.sum([MT['STELLARTYPEPSTMT2'] == 11 ])}") 
print(f"Number of secondary ONeWD (pst MT) :  {np.sum([MT['STELLARTYPEPSTMT2'] == 12 ])}")

print(f"Number of BH-MS binaries (pre MT): {np.sum(MASKPREMTBHMS1)}") 
print(f"Number of MS-BH binaries (pre MT): {np.sum(MASKPREMTBHMS2)}")
print(f"Number of BH-MS binaries (pst MT):  {np.sum(MASKPSTMTBHMS1)}") 
print(f"Number of MS-BH binaries (pst MT):  {np.sum(MASKPSTMTBHMS2)}") 
print(f"Number of BH-MS binaries (pst last MT):  {np.sum(MASKPSTMTBHMS1*MASKLASTMT)}") 
print(f"Number of MS-BH binaries (pst last MT):  {np.sum(MASKPSTMTBHMS2*MASKLASTMT)}") 
print(f"Number of systems with semimajor axis < 3.0 AU BH-MS binaries (pre MT):  {np.sum(MASKPREMTSEMAJ*MASKPREMTBHMS1)}") 
print(f"Number of systems with semimajor axis < 3.0 AU MS-BH binaries (pre MT):  {np.sum(MASKPREMTSEMAJ*MASKPREMTBHMS2)}") 
print(f"Number of systems with semimajor axis < 3.0 AU BH-MS binaries (pst MT):  {np.sum(MASKPSTMTSEMAJ*MASKPSTMTBHMS1)}") 
print(f"Number of systems with semimajor axis < 3.0 AU MS-BH binaries (pst MT):  {np.sum(MASKPSTMTSEMAJ*MASKPSTMTBHMS2)}") 
print(f"Number of systems with semimajor axis < 3.0 AU BH-MS binaries (pst last MT):  {np.sum(MASKPSTMTSEMAJ*MASKPSTMTBHMS1*MASKLASTMT)}") 
print(f"Number of systems with semimajor axis < 3.0 AU MS-BH binaries (pst last MT):  {np.sum(MASKPSTMTSEMAJ*MASKPSTMTBHMS2*MASKLASTMT)}")
print(f"Number of systems with orbital period < 70 days BH-MS binaries (pre MT):  {np.sum(MASKPREMTORBPER*MASKPREMTBHMS1)}")
print(f"Number of systems with orbital period < 70 days MS-BH binaries (pre MT):  {np.sum(MASKPREMTORBPER*MASKPREMTBHMS2)}")
print(f"Number of systems with orbital period < 70 days BH-MS binaries (pst MT):  {np.sum(MASKPSTMTORBPER*MASKPSTMTBHMS1)}")
print(f"Number of systems with orbital period < 70 days MS-BH binaries (pst MT):  {np.sum(MASKPSTMTORBPER*MASKPSTMTBHMS2)}")
print(f"Number of systems with orbital period < 70 days BH-MS binaries (pst last MT):  {np.sum(MASKPSTMTORBPER*MASKPSTMTBHMS1*MASKLASTMT)}")
print(f"Number of systems with orbital period < 70 days MS-BH binaries (pst last MT):  {np.sum(MASKPSTMTORBPER*MASKPSTMTBHMS2*MASKLASTMT)}")
print(f"Number of systems with orbital period < 30 days BH-MS binaries (pst MT):  {np.sum(MASKPSTMTORBPER_masuda*MASKPSTMTBHMS1)}")
print(f"Number of systems with orbital period < 30 days MS-BH binaries (pst MT):  {np.sum(MASKPSTMTORBPER_masuda*MASKPSTMTBHMS2)}")
print(f"Number of systems with orbital period < 30 days BH-MS binaries (pst last MT):  {np.sum(MASKPSTMTORBPER_masuda*MASKPSTMTBHMS1*MASKLASTMT)}")
print(f"Number of systems with orbital period < 30 days MS-BH binaries (pst last MT):  {np.sum(MASKPSTMTORBPER_masuda*MASKPSTMTBHMS2*MASKLASTMT)}")
print(f"Number of systems with semimajor axis < 0.4 AU BH-MS binaries (pst MT):  {np.sum(MASKPSTMTSEMAJ_masuda*MASKPSTMTBHMS1)}")
print(f"Number of systems with semimajor axis < 0.4 AU MS-BH binaries (pst MT):  {np.sum(MASKPSTMTSEMAJ_masuda*MASKPSTMTBHMS2)}")
print(f"Number of systems with semimajor axis < 0.4 AU BH-MS binaries (pst last MT):  {np.sum(MASKPSTMTSEMAJ_masuda*MASKPSTMTBHMS1*MASKLASTMT)}")
print(f"Number of systems with semimajor axis < 0.4 AU MS-BH binaries (pst last MT):  {np.sum(MASKPSTMTSEMAJ_masuda*MASKPSTMTBHMS2*MASKLASTMT)}")

print(f"Number of systems with orbital inclination (COSI) <  np.cos(np.deg2rad(89.99)) BH-MS binaries (pre MT):  {np.sum(MASKPREMTSEARCHABILITY*MASKPREMTBHMS1)}")
print(f"Number of systems with orbital inclination (COSI) <  np.cos(np.deg2rad(89.99)) MS-BH binaries (pre MT):  {np.sum(MASKPREMTSEARCHABILITY*MASKPREMTBHMS2)}")
print(f"Number of systems with orbital inclination (COSI) <  np.cos(np.deg2rad(89.99)) BH-MS binaries (pst MT):  {np.sum(MASKPSTMTSEARCHABILITY*MASKPSTMTBHMS1)}")
print(f"Number of systems with orbital inclination (COSI) <  np.cos(np.deg2rad(89.99)) MS-BH binaries (pst MT):  {np.sum(MASKPSTMTSEARCHABILITY*MASKPSTMTBHMS2)}")
print(f"Number of systems with orbital inclination (COSI) <  np.cos(np.deg2rad(89.99)) BH-MS binaries (pst last MT):  {np.sum(MASKPSTMTSEARCHABILITY*MASKPSTMTBHMS1*MASKLASTMT)}")
print(f"Number of systems with orbital inclination (COSI) <  np.cos(np.deg2rad(89.99)) MS-BH binaries (pst last MT):  {np.sum(MASKPSTMTSEARCHABILITY*MASKPSTMTBHMS2*MASKLASTMT)}")

SP_mask_hdu = fits.BinTableHDU(Table(data=[MASKSPBH1,MASKSPBH2, MASKSPunb,MASKSPdco,MASKSPmrgr,MASKMTinSP, ORBITALPERIODSP], 
                                    names=["MASKSPBH1","MASKSPBH2", "MASKSPunb","MASKSPdco","MASKSPmrgr","MASKMTinSP", "ORBITALPERIODSP"]))
MT_mask_hdu = fits.BinTableHDU(Table(data=[MASKPREMTBH1,MASKPREMTBH2,MASKPREMTBHNS1,MASKPREMTBHNS2,MASKPREMTBHMS1,MASKPREMTBHMS2,
                                            MASKPREMTnonBH, MASKPREMTSEMAJ, MASKPREMTnegCE, MASKPREMTposCE,MASKFIRSTMT,
                                            MASKPSTMTBH1,MASKPSTMTBH2,MASKPSTMTBHNS1,MASKPSTMTBHNS2,MASKPSTMTBHMS1,MASKPSTMTBHMS2,
                                            MASKPSTMTnonBH, MASKPSTMTSEMAJ,
                                            MASKPSTMTnegCE, MASKPSTMTposCE, MASKLASTMT, MASKTYPECHANGE1, MASKTYPECHANGE2,
                                            ORBITALPERIODPREMT,ORBITALPERIODPSTMT, MASKPREMTORBPER, COSIPREMT,MASKPSTMTORBPER, COSIPSTMT,
                                            MASKPREMTSEARCHABILITY, MASKPSTMTSEARCHABILITY, CEAFTERMT, MASKPREMTSEMAJ_masuda, MASKPSTMTSEMAJ_masuda,
                                            MASKPREMTORBPER_masuda, MASKPSTMTORBPER_masuda], 
                                    names=["MASKPREMTBH1","MASKPREMTBH2","MASKPREMTBHNS1","MASKPREMTBHNS2","MASKPREMTBHMS1","MASKPREMTBHMS2",
                                            "MASKPREMTnonBH", "MASKPREMTSEMAJ", "MASKPREMTnegCE", "MASKPREMTposCE","MASKFIRSTMT",
                                            "MASKPSTMTBH1","MASKPSTMTBH2","MASKPSTMTBHNS1","MASKPSTMTBHNS2", "MASKPSTMTBHMS1","MASKPSTMTBHMS2", 
                                            "MASKPSTMTnonBH", "MASKPSTMTSEMAJ",
                                            "MASKPSTMTnegCE", "MASKPSTMTposCE", "MASKLASTMT", "MASKTYPECHANGE1", "MASKTYPECHANGE2",
                                            "ORBITALPERIODPREMT", "ORBITALPERIODPSTMT","MASKPREMTORBPER", "COSIPREMT", "MASKPSTMTORBPER", "COSIPSTMT",
                                            "MASKPREMTSEARCHABILITY", "MASKPSTMTSEARCHABILITY","CEAFTERMT", "MASKPREMTSEMAJ_masuda", "MASKPSTMTSEMAJ_masuda",
                                            "MASKPREMTORBPER_masuda", "MASKPSTMTORBPER_masuda"]))
hdu.append(SP_mask_hdu)
hdu.append(MT_mask_hdu)
hdu.close()


e = datetime.now()
# Displays Time
current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
duration = e - s
print(" \n Finish time :", current_time)