####### PARENT SCRIPT: tertius_plotting_020525
#------> this version plots for contrained semimajor range < 0.5 AU for MTs
# Companion mass order issue solved from the previous version
# Includes binary data copnstruction to be used by the capture scripts
#%%
import os, sys
import numpy as np               # for handling arrays
import seaborn as sns
import numpy.ma as ma
import pandas as pd
import time                      # for finding computation time
from datetime import datetime
import matplotlib.pyplot as plt  #for plotting
import matplotlib
import seaborn as sns
from astropy.io import fits
from luxetenebrae import calculations as calc
from matplotlib.ticker import FuncFormatter, MultipleLocator, ScalarFormatter,LogFormatter, LogLocator
from luxetenebrae import utils as utils      # functions from utils.py
import datetime as dt
import logging            # for reading the COMPAS data
import time as t

sns.set_theme(style="ticks")
now = dt.datetime.now()
# Choose the modee to process
mode = 'WD_Enabled'

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


# Keeps corresponding numerical values 
SPs = []
MTs = []
CEs = []

STELLARTYPEZAMSSP1   =  []
STELLARTYPEZAMSSP2   =  []
STELLARTYPESP1   =  []
STELLARTYPESP2   =  []       
STELLARTYPEMT1   =  []
STELLARTYPEMT2   =  [] 

MASSZAMSSP1 = []
MASSZAMSSP2 = []
MASSMT1 = []
MASSMT2 = []

SEMIMAJORAXISZAMSSP = []
SEMIMAJORAXISMT = []

# Boolean values for masking
MASKSPBH1 = np.array([], dtype='bool')
MASKSPBH2 = np.array([], dtype='bool')
MASKMTBH1 = np.array([], dtype='bool')
MASKMTBH2 = np.array([], dtype='bool')

MASKSPunb = np.array([], dtype='bool')
MASKSPMTo = np.array([], dtype='bool')
MASKSPmrgr = np.array([], dtype='bool')

MASKMTinSP = np.array([], dtype='bool')
MASKMTBHNS1 = np.array([], dtype='bool')
MASKMTBHNS2 = np.array([], dtype='bool')

MASKMTnonBH =  np.array([], dtype='bool')
MASKMTSEMAJ = np.array([], dtype='bool')

# Dataframes to use with sns plotting package
data_SP = {}
SPdf = pd.DataFrame(data_SP)
data_MT = {}
MTdf = pd.DataFrame(data_MT)

# Reads files produced by ... and contructs the dataframes. SPs and MTs are separate since they have different sizes.
arrays_to_save = [SPs,MTs,CEs,STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,STELLARTYPEMT1,STELLARTYPEMT2,
                MASSZAMSSP1,MASSZAMSSP2,MASSMT1,MASSMT2,SEMIMAJORAXISZAMSSP,SEMIMAJORAXISMT,MASKSPBH1,MASKSPBH2,MASKMTBH1,MASKMTBH2,
                MASKSPunb,MASKSPMTo,MASKSPmrgr,MASKMTinSP,MASKMTBHNS1,MASKMTBHNS2, MASKMTnonBH, MASKMTSEMAJ]
filenames = ["SPs","MTs","CEs","STELLARTYPEZAMSSP1","STELLARTYPEZAMSSP2","STELLARTYPESP1","STELLARTYPESP2","STELLARTYPEMT1","STELLARTYPEMT2",
            "MASSZAMSSP1","MASSZAMSSP2","MASSMT1","MASSMT2","SEMIMAJORAXISZAMSSP","SEMIMAJORAXISMT","MASKSPBH1","MASKSPBH2","MASKMTBH1","MASKMTBH2",
            "MASKSPunb","MASKSPMTo","MASKSPmrgr","MASKMTinSP","MASKMTBHNS1","MASKMTBHNS2", "MASKMTnonBH", "MASKMTSEMAJ"]

fits_data = pathToFiles + "/stage_02_output.fits"
with fits.open(fits_data) as hdul:
    hdul.info()
    SP = hdul[1].data
    MT = hdul[2].data
    CE = hdul[3].data

fits_data = pathToFiles + "/stage_03_output.fits"
with fits.open(fits_data) as hdul:
    hdul.info()
    SP_mask = hdul[1].data
    MT_mask = hdul[2].data


events = MT['EVENTSMT']
print(min(events), max(events), np.mean(events))


systemSize = len(SP['SPs'])
print(f"System size: {systemSize}")
MTdf_pre_singleBH = pd.DataFrame()
MTdf_pre_tot = pd.DataFrame()
#Systems which primary is BH
MTdf_pre_singleBH["BH1"] = MT["MASSPREMT1"]*MT_mask["MASKPREMTBH1"]
MTdf_pre_singleBH["CP2"] = MT["MASSPREMT2"]*MT_mask["MASKPREMTBH1"]
MTdf_pre_singleBH["SA1"] = MT["SEMIMAJORAXISPREMT"]*MT_mask["MASKPREMTBH1"]
MTdf_pre_singleBH["CPR2"] = MT["RADIUSPREMT2"]*MT_mask["MASKPREMTBH1"]


#Systems which secondary is BH
MTdf_pre_singleBH["BH2"] = MT["MASSPREMT2"]*MT_mask["MASKPREMTBH2"]
MTdf_pre_singleBH["CP1"] = MT["MASSPREMT1"]*MT_mask["MASKPREMTBH2"]
MTdf_pre_singleBH["SA2"] = MT["SEMIMAJORAXISPREMT"]*MT_mask["MASKPREMTBH2"]
MTdf_pre_singleBH["CPR1"] = MT["RADIUSPREMT1"]*MT_mask["MASKPREMTBH2"]

########## UPDATE 7.11.24
#disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
MTdf_pre_tot["BlackHole"] = np.concatenate([MTdf_pre_singleBH["BH1"], MTdf_pre_singleBH["BH2"]])
MTdf_pre_tot["Companion"] = np.concatenate([MTdf_pre_singleBH["CP2"], MTdf_pre_singleBH["CP1"]])
MTdf_pre_tot["Semax"] = np.concatenate([MTdf_pre_singleBH["SA1"] , MTdf_pre_singleBH["SA2"] ])
MTdf_pre_tot["CPR"] = np.concatenate([MTdf_pre_singleBH["CPR2"] , MTdf_pre_singleBH["CPR1"]])

MTdf_pst_singleBH = pd.DataFrame()
MTdf_pst_tot = pd.DataFrame()
#Systems which primary is BH
MTdf_pst_singleBH["BH1"] = MT["MASSPSTMT1"]*MT_mask["MASKPSTMTBH1"]
MTdf_pst_singleBH["CP2"] = MT["MASSPSTMT2"]*MT_mask["MASKPSTMTBH1"]
MTdf_pst_singleBH["SA1"] = MT["SEMIMAJORAXISPSTMT"]*MT_mask["MASKPSTMTBH1"]
MTdf_pst_singleBH["OP1"] = MT_mask["ORBITALPERIODPSTMT"]*MT_mask["MASKPSTMTBH1"]
MTdf_pst_singleBH["CPR2"] = MT["RADIUSPSTMT2"]*MT_mask["MASKPSTMTBH1"]

#Systems which secondary is BH
MTdf_pst_singleBH["BH2"] = MT["MASSPSTMT2"]*MT_mask["MASKPSTMTBH2"]
MTdf_pst_singleBH["CP1"] = MT["MASSPSTMT1"]*MT_mask["MASKPSTMTBH2"]
MTdf_pst_singleBH["SA2"] = MT["SEMIMAJORAXISPSTMT"]*MT_mask["MASKPSTMTBH2"]
MTdf_pst_singleBH["OP2"] = MT_mask["ORBITALPERIODPSTMT"]*MT_mask["MASKPSTMTBH2"]
MTdf_pst_singleBH["CPR1"] = MT["RADIUSPSTMT1"]*MT_mask["MASKPSTMTBH2"]
########## UPDATE 7.11.24
#disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
MTdf_pst_tot["BlackHole"] = np.concatenate([MTdf_pst_singleBH["BH1"], MTdf_pst_singleBH["BH2"]])
MTdf_pst_tot["Companion"] = np.concatenate([MTdf_pst_singleBH["CP2"], MTdf_pst_singleBH["CP1"]])
MTdf_pst_tot["Semax"] = np.concatenate([MTdf_pst_singleBH["SA1"] , MTdf_pst_singleBH["SA2"] ])
MTdf_pst_tot["OP"] = np.concatenate([MTdf_pst_singleBH["OP1"] , MTdf_pst_singleBH["OP2"] ])
MTdf_pst_tot["CPR"] = np.concatenate([MTdf_pst_singleBH["CPR2"] , MTdf_pst_singleBH["CPR1"]])

print('Companion: ', np.min(MTdf_pst_tot["Companion"]), np.max(MTdf_pst_tot["Companion"]))

change = (MT_mask['MASKTYPECHANGE1'] == 1) | (MT_mask['MASKTYPECHANGE2'] == 1)
MTdf_pst_tot["BlackHole_changemsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*change, MTdf_pst_singleBH["BH2"]*change])
MTdf_pst_tot["Companion_changemsk"] = np.concatenate([MTdf_pst_singleBH["CP2"]*change, MTdf_pst_singleBH["CP1"]*change])
MTdf_pst_tot["Semax_changemsk"] = np.concatenate([MTdf_pst_singleBH["SA1"]*change , MTdf_pst_singleBH["SA2"]*change ])  

first = MT_mask['MASKFIRSTMT']
MTdf_pre_tot["BlackHole_firstmsk"] = np.concatenate([MTdf_pre_singleBH["BH1"]*first, MTdf_pre_singleBH["BH2"]*first])
MTdf_pre_tot["Companion_firstmsk"] = np.concatenate([MTdf_pre_singleBH["CP2"]*first, MTdf_pre_singleBH["CP1"]*first])
MTdf_pre_tot["Semax_firstmsk"] = np.concatenate([MTdf_pre_singleBH["SA1"]*first , MTdf_pre_singleBH["SA2"]*first ]) 

last = MT_mask['MASKLASTMT']
MTdf_pst_tot["BlackHole_lastmsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*last, MTdf_pst_singleBH["BH2"]*last])
MTdf_pst_tot["Companion_lastmsk"] = np.concatenate([MTdf_pst_singleBH["CP2"]*last, MTdf_pst_singleBH["CP1"]*last])
MTdf_pst_tot["Semax_lastmsk"] = np.concatenate([MTdf_pst_singleBH["SA1"]*last , MTdf_pst_singleBH["SA2"]*last ])

print("Checking for NaN values in MTdf_pst_tot['Semax_lastmsk']:", np.isnan(*calc.mylog(MTdf_pst_tot["Semax_lastmsk"],MTdf_pst_tot["BlackHole_lastmsk"])).any())


mainsequence_pre = (MT_mask['MASKPREMTBHMS1'] | MT_mask['MASKPREMTBHMS2']) & first
MTdf_pre_tot["BlackHole_mainpremsk"] = np.concatenate([MTdf_pre_singleBH["BH1"]*mainsequence_pre, MTdf_pre_singleBH["BH2"]*mainsequence_pre])
MTdf_pre_tot["Companion_mainpremsk"] = np.concatenate([MTdf_pre_singleBH["CP2"]*mainsequence_pre, MTdf_pre_singleBH["CP1"]*mainsequence_pre])
MTdf_pre_tot["Semax_mainpremsk"] = np.concatenate([MTdf_pre_singleBH["SA1"]*mainsequence_pre , MTdf_pre_singleBH["SA2"]*mainsequence_pre ])

mainsequence_pst = (MT_mask['MASKPSTMTBHMS1'] | MT_mask['MASKPSTMTBHMS2']) & last
MTdf_pst_tot["BlackHole_mainpstmsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*mainsequence_pst, MTdf_pst_singleBH["BH2"]*mainsequence_pst])
MTdf_pst_tot["Companion_mainpstmsk"] = np.concatenate([MTdf_pst_singleBH["CP2"]*mainsequence_pst, MTdf_pst_singleBH["CP1"]*mainsequence_pst])
MTdf_pst_tot["Semax_mainpstmsk"] = np.concatenate([MTdf_pst_singleBH["SA1"]*mainsequence_pst , MTdf_pst_singleBH["SA2"]*mainsequence_pst ])
MTdf_pst_tot["OP_mainpstmsk"] = np.concatenate([MTdf_pst_singleBH["OP1"]*mainsequence_pst , MTdf_pst_singleBH["OP2"]*mainsequence_pst ])
MTdf_pst_tot["CompanionR_mainpstmsk"] = np.concatenate([MTdf_pst_singleBH["CPR2"]*mainsequence_pst , MTdf_pst_singleBH["CPR1"]*mainsequence_pst ])
MTdf_pst_tot["BlackHoleR_mainpstmsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*mainsequence_pst, MTdf_pst_singleBH["BH2"]*mainsequence_pst])
mainsequence_pre_Tlim = mainsequence_pre & MT_mask['MASKPREMTORBPER']
MTdf_pre_tot["BlackHole_mainpreTlimmsk"] = np.concatenate([MTdf_pre_singleBH["BH1"]*mainsequence_pre_Tlim, MTdf_pre_singleBH["BH2"]*mainsequence_pre_Tlim])
MTdf_pre_tot["Companion_mainpreTlimmsk"] = np.concatenate([MTdf_pre_singleBH["CP2"]*mainsequence_pre_Tlim, MTdf_pre_singleBH["CP1"]*mainsequence_pre_Tlim])
MTdf_pre_tot["Semax_mainpreTlimmsk"] = np.concatenate([MTdf_pre_singleBH["SA1"]*mainsequence_pre_Tlim , MTdf_pre_singleBH["SA2"]*mainsequence_pre_Tlim ])

mainsequence_pst_Tlim = mainsequence_pst & MT_mask['MASKPSTMTORBPER']
MTdf_pst_tot["BlackHole_mainpstTlimmsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*mainsequence_pst_Tlim, MTdf_pst_singleBH["BH2"]*mainsequence_pst_Tlim])
MTdf_pst_tot["Companion_mainpstTlimmsk"] = np.concatenate([MTdf_pst_singleBH["CP2"]*mainsequence_pst_Tlim, MTdf_pst_singleBH["CP1"]*mainsequence_pst_Tlim])
MTdf_pst_tot["Semax_mainpstTlimmsk"] = np.concatenate([MTdf_pst_singleBH["SA1"]*mainsequence_pst_Tlim , MTdf_pst_singleBH["SA2"]*mainsequence_pst_Tlim ])
MTdf_pst_tot["OP_mainpstTlimmsk"] = np.concatenate([MTdf_pst_singleBH["OP1"]*mainsequence_pst_Tlim , MTdf_pst_singleBH["OP2"]*mainsequence_pst_Tlim ])

print('Companion: ', np.min(MTdf_pst_tot["Companion_mainpstTlimmsk"]), np.max(MTdf_pst_tot["Companion_mainpstTlimmsk"]))

searchability = calc.searchability(MTdf_pst_singleBH["SA1"])
N = np.sum(searchability)
searchability_index = N/len(searchability)
# searchability_index = np.sum(MT_mask['MASKPSTMTSEARCHABILITY'])/len(MT_mask['MASKPSTMTSEARCHABILITY'])
print('Searchability index:', searchability_index)
mainsequence_pst_alim_masuda = mainsequence_pst & MT_mask['MASKPSTMTSEMAJ_masuda'] 
MTdf_pst_tot["BlackHole_mainpstalimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*mainsequence_pst_alim_masuda, MTdf_pst_singleBH["BH2"]*mainsequence_pst_alim_masuda])
MTdf_pst_tot["Companion_mainpstalimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["CP2"]*mainsequence_pst_alim_masuda, MTdf_pst_singleBH["CP1"]*mainsequence_pst_alim_masuda])
MTdf_pst_tot["Semax_mainpstalimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["SA1"]*mainsequence_pst_alim_masuda , MTdf_pst_singleBH["SA2"]*mainsequence_pst_alim_masuda ])
MTdf_pst_tot["OP_mainpstalimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["OP1"]*mainsequence_pst_alim_masuda , MTdf_pst_singleBH["OP2"]*mainsequence_pst_alim_masuda ])

# print('Number of systems past last MT event with Orbital period less than 30 days:', np.sum(mainsequence_pst_alim_masuda))
print('Number of searchable systems past last MT event with semaj less than 0.4 au:', searchability_index*np.sum(mainsequence_pst_alim_masuda))
print('Number of systems past last with semaj less than 0.4 au:', np.sum(mainsequence_pst_alim_masuda))

mainsequence_pst_Tlim_masuda = mainsequence_pst & MT_mask['MASKPSTMTORBPER_masuda']
MTdf_pst_tot["BlackHole_mainpstTlimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*mainsequence_pst_Tlim_masuda, MTdf_pst_singleBH["BH2"]*mainsequence_pst_Tlim_masuda])
MTdf_pst_tot["Companion_mainpstTlimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["CP2"]*mainsequence_pst_Tlim_masuda, MTdf_pst_singleBH["CP1"]*mainsequence_pst_Tlim_masuda])
MTdf_pst_tot["Semax_mainpstTlimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["SA1"]*mainsequence_pst_Tlim_masuda , MTdf_pst_singleBH["SA2"]*mainsequence_pst_Tlim_masuda ])
MTdf_pst_tot["OP_mainpstTlimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["OP1"]*mainsequence_pst_Tlim_masuda , MTdf_pst_singleBH["OP2"]*mainsequence_pst_Tlim_masuda ])
MTdf_pst_tot["CPR_mainpstTlimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["CPR2"]*mainsequence_pst_Tlim_masuda , MTdf_pst_singleBH["CPR1"]*mainsequence_pst_Tlim_masuda ])

print('Number of searchable systems past last MT event with Orbital period less than 30 days:', searchability_index* np.sum(mainsequence_pst_Tlim_masuda))
print('Number of systems past last MT event with Orbital period less than 30 days:', np.sum(mainsequence_pst_Tlim_masuda))
print('Number of systems past last MT with Orbital period less than 30 days:', np.sum((MT_mask['MASKPSTMTORBPER_masuda']) * last))
print('Number of systems past last MT with orbital period less than 70 days:', np.sum((MT_mask['MASKPSTMTORBPER']) * last))
print('Number of bhms past last MT with Orbital period less than 30 days:', np.sum((MT_mask['MASKPSTMTORBPER_masuda']) * mainsequence_pst))
print('Number of bhms past last MT with orbital period less than 70 days:', np.sum((MT_mask['MASKPSTMTORBPER']) * mainsequence_pst))
print('max Orbital period in days masuda:', np.max(MTdf_pst_tot["OP_mainpstTlimmasudamsk"]))
print('mean Orbital period in days masuda:', np.mean(MTdf_pst_tot["OP_mainpstTlimmasudamsk"]))
print('max semimajor axis in au masuda:', np.max(MTdf_pst_tot["Semax_mainpstTlimmasudamsk"]))
print('max orbital period in days:', np.max(MTdf_pst_tot["OP_mainpstTlimmsk"]))
print('mean orbital period in days:', np.mean(MTdf_pst_tot["OP_mainpstTlimmsk"]))
print('max semimajor axis in au:', np.max(MTdf_pst_tot["Semax_mainpstTlimmsk"]))
print('min semimajor axis in au:', np.min(MTdf_pst_tot["Semax_mainpstTlimmsk"]), 'number of systems:', len(MTdf_pst_tot["Semax_mainpstTlimmsk"]))
print('min non zero semimajor axis in au:', np.min([x for x in MTdf_pst_tot["Semax_mainpstTlimmsk"] if x > 0.0]), 'number of systems:', len([x for x in MTdf_pst_tot["Semax_mainpstTlimmsk"] if x > 0.0]))
print('star radius and bh radius where the SA is minimum:', MTdf_pst_tot["CompanionR_mainpstmsk"][np.where(MTdf_pst_tot["Semax_mainpstTlimmsk"]== np.min([x for x in MTdf_pst_tot["Semax_mainpstTlimmsk"] if x > 0.0]))[0][0]] )

SPdf_singleBH = pd.DataFrame()
SPdf_tot = pd.DataFrame()
    #Systems which primary is BH
SPdf_singleBH["BH1"] = SP["MASSZAMSSP1"]*SP_mask["MASKSPBH1"]
SPdf_singleBH["CP2"] = SP["MASSZAMSSP2"]*SP_mask["MASKSPBH1"]
SPdf_singleBH["SA1"] = SP["SEMIMAJORAXISZAMSSP"]*SP_mask["MASKSPBH1"]

#Systems which secondary is BH
SPdf_singleBH["BH2"] = SP["MASSZAMSSP2"]*SP_mask["MASKSPBH2"]
SPdf_singleBH["CP1"] = SP["MASSZAMSSP1"]*SP_mask["MASKSPBH2"]
SPdf_singleBH["SA2"] = SP["SEMIMAJORAXISZAMSSP"]*SP_mask["MASKSPBH2"]

SPdf_singleMS = pd.DataFrame()
SPdf_totMS = pd.DataFrame()
#Sytems which primary is MS
SPdf_singleMS["MS1"] = SP["MASSZAMSSP1"]*SP_mask["MASKSPMS1"]
SPdf_singleMS["CP2"] = SP["MASSZAMSSP2"]*SP_mask["MASKSPMS1"]
SPdf_singleMS["SA1"] = SP["SEMIMAJORAXISZAMSSP"]*SP_mask["MASKSPMS1"]
SPdf_singleMS["MSR1"] = SP["RADIUSSP1"]*SP_mask["MASKSPMS1"]
SPdf_singleMS["CPR2"] = SP["RADIUSSP2"]*SP_mask["MASKSPMS1"]
SPdf_singleMS["OP1"] = SP_mask["ORBITALPERIODSP"]*SP_mask["MASKSPMS1"]
#Systems which secondary is MS
SPdf_singleMS["MS2"] = SP["MASSZAMSSP2"]*SP_mask["MASKSPMS2"]
SPdf_singleMS["CP1"] = SP["MASSZAMSSP1"]*SP_mask["MASKSPMS2"]
SPdf_singleMS["SA2"] = SP["SEMIMAJORAXISZAMSSP"]*SP_mask["MASKSPMS2"]
SPdf_singleMS["MSR2"] = SP["RADIUSSP2"]*SP_mask["MASKSPMS2"]
SPdf_singleMS["CPR1"] = SP["RADIUSSP1"]*SP_mask["MASKSPMS2"]
SPdf_singleMS["OP2"] = SP_mask["ORBITALPERIODSP"]*SP_mask["MASKSPMS2"]

########## UPDATE 7.11.24
#disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
SPdf_tot["BlackHole"] = np.concatenate([SPdf_singleBH["BH1"], SPdf_singleBH["BH2"]])
SPdf_tot["Companion"] = np.concatenate([SPdf_singleBH["CP2"], SPdf_singleBH["CP1"]])
SPdf_tot["Semax"] = np.concatenate([SPdf_singleBH["SA1"] , SPdf_singleBH["SA2"] ])

SPdf_totMS["MainSequence"] = np.concatenate([SPdf_singleMS["MS1"], SPdf_singleMS["MS2"]])
SPdf_totMS["Companion"] = np.concatenate([SPdf_singleMS["CP2"], SPdf_singleMS["CP1"]])
SPdf_totMS["Semax"] = np.concatenate([SPdf_singleMS["SA1"] , SPdf_singleMS["SA2"] ])
SPdf_totMS["OP"] = np.concatenate([SPdf_singleMS["OP1"] , SPdf_singleMS["OP2"] ])
SPdf_totMS["MainSequenceR"] = np.concatenate([SPdf_singleMS["MSR1"] , SPdf_singleMS["MSR2"] ])
SPdf_totMS["CompanionR"] = np.concatenate([SPdf_singleMS["CPR2"] , SPdf_singleMS["CPR1"] ])
##########

    #Systems which primary is BH (MTs in SP)
SPdf_singleBH["BH1_MT"] = SP["MASSZAMSSP1"]*SP_mask["MASKSPBH1"]*SP_mask["MASKMTinSP"]
SPdf_singleBH["CP2_MT"] = SP["MASSZAMSSP2"]*SP_mask["MASKSPBH1"]*SP_mask["MASKMTinSP"]
SPdf_singleBH["SA1_MT"] = SP["SEMIMAJORAXISZAMSSP"]*SP_mask["MASKSPBH1"]*SP_mask["MASKMTinSP"]

#Systems which secondary is BH (MTs in SP)
SPdf_singleBH["BH2_MT"] = SP["MASSZAMSSP2"]*SP_mask["MASKSPBH2"]*SP_mask["MASKMTinSP"]
SPdf_singleBH["CP1_MT"] = SP["MASSZAMSSP1"]*SP_mask["MASKSPBH2"]*SP_mask["MASKMTinSP"]
SPdf_singleBH["SA2_MT"] = SP["SEMIMAJORAXISZAMSSP"]*SP_mask["MASKSPBH2"]*SP_mask["MASKMTinSP"]

########## UPDATE 7.11.24
#disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary (MTs in SP)
SPdf_tot["BlackHole_MT"] = np.concatenate([SPdf_singleBH["BH1_MT"], SPdf_singleBH["BH2_MT"]])
SPdf_tot["Companion_MT"] = np.concatenate([SPdf_singleBH["CP2_MT"], SPdf_singleBH["CP1_MT"]])
SPdf_tot["Semax_MT"] = np.concatenate([SPdf_singleBH["SA1_MT"] , SPdf_singleBH["SA2_MT"] ])
##########



# Prepare the data
companion_mass = MT['MASSPSTMT2'] * MT_mask['MASKPSTMTMS1'] * MT_mask['MASKLASTMT']
msstar_mass = MT['MASSPSTMT1'] * MT_mask['MASKPSTMTMS1']* MT_mask['MASKLASTMT']
companion_radius = MT['RADIUSPSTMT2'] * MT_mask['MASKPSTMTMS1']* MT_mask['MASKLASTMT']
msstar_radius = MT['RADIUSPSTMT1'] * MT_mask['MASKPSTMTMS1']* MT_mask['MASKLASTMT']
semaj = MT['SEMIMAJORAXISPSTMT'] * MT_mask['MASKPSTMTMS1']* MT_mask['MASKLASTMT']
orbital_period = MT_mask['ORBITALPERIODPSTMT'] * MT_mask['MASKPSTMTMS1']* MT_mask['MASKLASTMT']
print(np.sum(MT_mask['MASKPSTMTMS1']* MT_mask['MASKLASTMT']))
# Create a FITS table
col1 = fits.Column(name="Companion_mass", array=companion_mass, format="E")
col2 = fits.Column(name="MSstar_mass", array=msstar_mass, format="E")
col3 = fits.Column(name="Companion_radius", array=companion_radius, format="E")
col4 = fits.Column(name="MSstar_radius", array=msstar_radius, format="E")
col5 = fits.Column(name="Semimajor_axis", array=semaj, format="E")
col6 = fits.Column(name="Orbital_period", array=orbital_period, format="E")
hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
# Save the FITS file
output_fits_path = "/data/a.saricaoglu/repo/COMPAS/Files/msstarvsx2_pstlastMT.fits"
hdu.writeto(output_fits_path, overwrite=True)
print(f"Data saved to {output_fits_path}")  

companion_mass = MT['MASSPREMT2'] * MT_mask['MASKPREMTMS1'] * MT_mask['MASKFIRSTMT']
msstar_mass = MT['MASSPREMT1'] * MT_mask['MASKPREMTMS1']* MT_mask['MASKFIRSTMT']
companion_radius = MT['RADIUSPREMT2'] * MT_mask['MASKPREMTMS1']* MT_mask['MASKFIRSTMT']
msstar_radius = MT['RADIUSPREMT1'] * MT_mask['MASKPREMTMS1']* MT_mask['MASKFIRSTMT']
semaj = MT['SEMIMAJORAXISPREMT'] * MT_mask['MASKPREMTMS1']* MT_mask['MASKFIRSTMT']
orbital_period = MT_mask['ORBITALPERIODPREMT'] * MT_mask['MASKPREMTMS1']* MT_mask['MASKFIRSTMT']
print(np.sum(MT_mask['MASKPREMTMS1']* MT_mask['MASKFIRSTMT']))
# Create a FITS table
col1 = fits.Column(name="Companion_mass", array=companion_mass, format="E")
col2 = fits.Column(name="MSstar_mass", array=msstar_mass, format="E")
col3 = fits.Column(name="Companion_radius", array=companion_radius, format="E")
col4 = fits.Column(name="MSstar_radius", array=msstar_radius, format="E")
col5 = fits.Column(name="Semimajor_axis", array=semaj, format="E")
col6 = fits.Column(name="Orbital_period", array=orbital_period, format="E")
hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
# Save the FITS file
output_fits_path = "/data/a.saricaoglu/repo/COMPAS/Files/msstarvsx_prefirstMT.fits"
hdu.writeto(output_fits_path, overwrite=True)
print(f"Data saved to {output_fits_path}")    



# # Prepare the data
# companion_mass = MTdf_pst_tot["Companion_mainpstmsk"]
# blackhole_mass = MTdf_pst_tot["BlackHole_mainpstmsk"]
# companion_radius = MTdf_pst_tot["CompanionR_mainpstmsk"]
# blackhole_radius = MTdf_pst_tot["BlackHoleR_mainpstmsk"]
# semaj = MTdf_pst_tot["Semax_mainpstmsk"]
# orbital_period = MTdf_pst_tot["OP_mainpstmsk"]

# # Create a FITS table
# col1 = fits.Column(name="Companion_mass", array=companion_mass, format="E")
# col2 = fits.Column(name="Blackhole_mass", array=blackhole_mass, format="E")
# col3 = fits.Column(name="Companion_radius", array=companion_radius, format="E")
# col4 = fits.Column(name="Blackhole_radius", array=blackhole_radius, format="E")
# col5 = fits.Column(name="Semimajor_axis", array=semaj, format="E")
# col6 = fits.Column(name="Orbital_period", array=orbital_period, format="E")
# hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])

# # Save the FITS file
# output_fits_path = "/data/a.saricaoglu/repo/COMPAS/Files/bhvscomp.fits"
# hdu.writeto(output_fits_path, overwrite=True)
# print(f"Data saved to {output_fits_path}")

# # Prepare the data
# companion_mass = SPdf_totMS["Companion"]
# msstar_mass = SPdf_totMS["MainSequence"]
# companion_radius = SPdf_totMS["CompanionR"]
# msstar_radius = SPdf_totMS["MainSequenceR"]
# semaj = SPdf_totMS["Semax"]
# orbital_period = SPdf_totMS["OP"]
# # Create a FITS table
# col1 = fits.Column(name="Companion_mass", array=companion_mass, format="E")
# col2 = fits.Column(name="MSstar_mass", array=msstar_mass, format="E")
# col3 = fits.Column(name="Companion_radius", array=companion_radius, format="E")
# col4 = fits.Column(name="MSstar_radius", array=msstar_radius, format="E")
# col5 = fits.Column(name="Semimajor_axis", array=semaj, format="E")
# col6 = fits.Column(name="Orbital_period", array=orbital_period, format="E")
# hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
# # Save the FITS file
# print(f"length of msstar:", len(msstar_mass))
# output_fits_path = "/data/a.saricaoglu/repo/COMPAS/Files/msstarvscomp_SP.fits"
# hdu.writeto(output_fits_path, overwrite=True)
# print(f"Data saved to {output_fits_path}")




