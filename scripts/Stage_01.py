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
from astropy.io import fits
from astropy.table import Table
from astropy import constants as const
from astropy import units as u
from luxetenebrae import utils as utils


now = dt.datetime.now()
# Choose the mode to process
mode = 'WD_Enabled'

# Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print("COMPAS_ROOT_DIR:", compasRootDir)
print(sys.path)

pathToData = os.path.join(compasRootDir, "luxetenebrae/runs", mode)
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
runs= [x[0] for x in os.walk(pathToData) if "COMPAS" in x[0]]
data_outputs = []
for run in runs:
    out = [f for f in os.listdir(run) if ".h5" in f]
    try:
        data = h5.File(run + "/" + out[0])
        print("Reading file: ", run + "/" + out[0])
        data_outputs.append(data)
    except:
        continue



fit_filename = pathToFiles + "/stage_01_output.fits"
hdu_pr = fits.PrimaryHDU()
hdu_pr.writeto(fit_filename, overwrite=True)
hdu = fits.open(fit_filename, mode='update')

SPs = []

EVENTSMT = []
EVOLUTIONSTATSP = []

STELLARTYPEZAMSSP1   =  []
STELLARTYPEZAMSSP2   =  []
STELLARTYPESP1   =  []
STELLARTYPESP2   =  []   
MASSZAMSSP1 = []
MASSZAMSSP2 = []
EVOLUTIONSTATSP = []
SEMIMAJORAXISZAMSSP = []
ECCENTRICITYZAMSSP = []
CIRCULARIZATIONTIME = []
RADIUSSP1 = []
RADIUSSP2 = []

MTs = []
EVENTSMT = []
STELLARTYPEPREMT1   =  []
STELLARTYPEPREMT2   =  [] 
STELLARTYPEPSTMT1   =  []
STELLARTYPEPSTMT2   =  [] 
MASSPREMT1 = []
MASSPREMT2 = []
MASSPSTMT1 = []
MASSPSTMT2 = []

SEMIMAJORAXISZAMSSP = []
SEMIMAJORAXISPREMT = []
SEMIMAJORAXISPSTMT = []
ECCENTRICITYPREMT = []
ECCENTRICITYPSTMT = []
RADIUSPREMT1 = []
RADIUSPREMT2 = []
RADIUSPSTMT1 = []
RADIUSPSTMT2 = []
TIMEPREMT = []
TIMEPSTMT = []
CEAFTERMT = []
MASSTRANSFERHISTORY = []

CEs = []
EVENTSCE = []
STELLARTYPEPRECE1   =  []
STELLARTYPEPRECE2   =  [] 
STELLARTYPEPSTCE1   =  []
STELLARTYPEPSTCE2   =  [] 
MASSPRECE1 = []
MASSPRECE2 = []
MASSPSTCE1 = []
MASSPSTCE2 = []
SEMIMAJORAXISPRECE = []
SEMIMAJORAXISPSTCE = []
RADIUSPRECE1 = []
RADIUSPRECE2 = []
RADIUSPSTCE1 = []
RADIUSPSTCE2 = []
COSIPRECE = []
COSIPSTCE = []
TIMECE = []

print(f'Data outputs :, {data_outputs}')

i = 0

for Data in data_outputs:

    c = dt.datetime.now()
    current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')

    SP = Data['BSE_System_Parameters']
    # print(SP.keys())



    seedsSP = SP['SEED'][()]
    statusSP = SP['Evolution_Status'][()]

    print(f"Batch (of {len(seedsSP)} sys.)" + str(i) +  " start time :", current_time)

    stellarTypeZamsSP1   =  SP['Stellar_Type@ZAMS(1)'][()]
    stellarTypeZamsSP2   =  SP['Stellar_Type@ZAMS(2)'][()]
    stellarTypeSP1   =  SP['Stellar_Type(1)'][()]
    stellarTypeSP2   =  SP['Stellar_Type(2)'][()]  

    massZamsSP1 = SP['Mass@ZAMS(1)'][()] 
    massZamsSP2 = SP['Mass@ZAMS(2)'][()]

    semimajorAxisZamsSP = SP['SemiMajorAxis@ZAMS'][()] # in AU
    eccentricityZamsSP = SP['Eccentricity@ZAMS'][()]

    radiusSP1 = massZamsSP1**(0.8) # in R_sun
    radiusSP2 = massZamsSP2**(0.8) # in R_sun
    print(radiusSP1)


    

    SPs.extend(seedsSP)
    EVOLUTIONSTATSP.extend(statusSP)

    STELLARTYPEZAMSSP1.extend(stellarTypeZamsSP1)
    STELLARTYPEZAMSSP2.extend(stellarTypeZamsSP2)
    STELLARTYPESP1.extend(stellarTypeSP1)
    STELLARTYPESP2.extend(stellarTypeSP2)   

    MASSZAMSSP1.extend(massZamsSP1)
    MASSZAMSSP2.extend(massZamsSP2)

    SEMIMAJORAXISZAMSSP.extend(semimajorAxisZamsSP)
    ECCENTRICITYZAMSSP.extend(eccentricityZamsSP)

    RADIUSSP1.extend(radiusSP1)
    RADIUSSP2.extend(radiusSP2)



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
    
    massTransferhistory = MT['CEE>MT'][()]

    radiuspreMT1 = MT['Radius(1)<MT'][()]
    radiuspstMT1 = MT['Radius(1)>MT'][()] 
    radiuspreMT2 = MT['Radius(2)<MT'][()]
    radiuspstMT2 = MT['Radius(2)>MT'][()]

    

    CEafterMT = MT['CEE>MT'][()]

    MTs.extend(seedsMT)
    EVENTSMT.extend(eventsMT)
    
    STELLARTYPEPREMT1.extend(stellarTypepreMT1)
    STELLARTYPEPREMT2.extend(stellarTypepreMT2) 
    STELLARTYPEPSTMT1.extend(stellarTypepstMT1)
    STELLARTYPEPSTMT2.extend(stellarTypepstMT2) 

    MASSPREMT1.extend(masspreMT1)
    MASSPREMT2.extend(masspreMT2)
    MASSPSTMT1.extend(masspstMT1)
    MASSPSTMT2.extend(masspstMT2)

    SEMIMAJORAXISPREMT.extend(semimajorAxispreMT)
    SEMIMAJORAXISPSTMT.extend(semimajorAxispstMT)

    ECCENTRICITYPREMT.extend(eccentricitypreMT)
    ECCENTRICITYPSTMT.extend(eccentricitypstMT)

    TIMEPREMT.extend(timepreMT)
    TIMEPSTMT.extend(timepstMT)

    MASSTRANSFERHISTORY.extend(massTransferhistory)
    


    CEAFTERMT.extend(CEafterMT)

    RADIUSPREMT1.extend(radiuspreMT1)
    RADIUSPREMT2.extend(radiuspreMT2)
    RADIUSPSTMT1.extend(radiuspstMT1)
    RADIUSPSTMT2.extend(radiuspstMT2)



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
    circulartizationTime = CE['Tau_Circ'][()]
    
    radiuspreCE1 = CE['Radius(1)<CE'][()] # in R_sun
    radiuspstCE1 = CE['Radius(1)>CE'][()] 
    radiuspreCE2 = CE['Radius(2)<CE'][()]
    radiuspstCE2 = CE['Radius(2)>CE'][()] 

    CEs.extend(seedsCE)
    EVENTSCE.extend(eventsCE)

    STELLARTYPEPRECE1.extend(stellarTypepreCE1)
    STELLARTYPEPRECE2.extend(stellarTypepreCE2) 
    STELLARTYPEPSTCE1.extend(stellarTypepstCE1)
    STELLARTYPEPSTCE2.extend(stellarTypepstCE2) 

    MASSPRECE1.extend(masspreCE1)
    MASSPRECE2.extend(masspreCE2)
    MASSPSTCE1.extend(masspstCE1)
    MASSPSTCE2.extend(masspstCE2)

    SEMIMAJORAXISPRECE.extend(semimajorAxispreCE)
    SEMIMAJORAXISPSTCE.extend(semimajorAxispstCE)

    CIRCULARIZATIONTIME.extend(circulartizationTime)


    # PreMT refers to the mass just before the first MT event, while pstMT refers to the mass just after the last MT event.


    c = dt.datetime.now()
    current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
    print("Batch " + str(i) +  " end time :", current_time)
    i = i+1
Data.close()

arrays_to_save = [ ]
filenames = []
checklist = [MTs,STELLARTYPEPREMT1,STELLARTYPEPREMT2,STELLARTYPEPSTMT1,STELLARTYPEPSTMT2,
                                    MASSPREMT1,MASSPREMT2,MASSPSTMT1,MASSPSTMT2,SEMIMAJORAXISPREMT, SEMIMAJORAXISPSTMT,
                                    EVENTSMT,MASSTRANSFERHISTORY,RADIUSPREMT1, RADIUSPREMT2, RADIUSPSTMT1, RADIUSPSTMT2,
                                    TIMEPREMT, TIMEPSTMT, CEAFTERMT]
for e in checklist:
    print(len(e))

SP_hdu = fits.BinTableHDU(Table(data=[SPs, STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,
                                    MASSZAMSSP1,MASSZAMSSP2,SEMIMAJORAXISZAMSSP, EVOLUTIONSTATSP, ECCENTRICITYZAMSSP, RADIUSSP1, RADIUSSP2], 
                                names=["SPs","STELLARTYPEZAMSSP1", "STELLARTYPEZAMSSP2", "STELLARTYPESP1", "STELLARTYPESP2",
                                    "MASSZAMSSP1", "MASSZAMSSP2",  "SEMIMAJORAXISZAMSSP", "EVOLUTIONSTATSP", "ECCENTRICITYZAMSSP", "RADIUSSP1", "RADIUSSP2"]))

MT_hdu = fits.BinTableHDU(Table(data=[MTs,STELLARTYPEPREMT1,STELLARTYPEPREMT2,STELLARTYPEPSTMT1,STELLARTYPEPSTMT2,
                                    MASSPREMT1,MASSPREMT2,MASSPSTMT1,MASSPSTMT2,SEMIMAJORAXISPREMT, SEMIMAJORAXISPSTMT,
                                    EVENTSMT,MASSTRANSFERHISTORY,RADIUSPREMT1, RADIUSPREMT2, RADIUSPSTMT1, RADIUSPSTMT2,
                                    TIMEPREMT, TIMEPSTMT, CEAFTERMT, ECCENTRICITYPREMT, ECCENTRICITYPSTMT], 
                                names=["MTs","STELLARTYPEPREMT1", "STELLARTYPEPREMT2", "STELLARTYPEPSTMT1", "STELLARTYPEPSTMT2",
                                    "MASSPREMT1", "MASSPREMT2", "MASSPSTMT1", "MASSPSTMT2", "SEMIMAJORAXISPREMT", "SEMIMAJORAXISPSTMT",
                                    "EVENTSMT","MASSTRANSFERHISTORY", "RADIUSPREMT1", "RADIUSPREMT2", "RADIUSPSTMT1", "RADIUSPSTMT2",
                                    "TIMEPREMT", "TIMEPSTMT","CEAFTERMT", "ECCENTRICITYPREMT", "ECCENTRICITYPSTMT"]))
CE_hdu = fits.BinTableHDU(Table(data=[CEs, EVENTSCE, STELLARTYPEPRECE1, STELLARTYPEPRECE2, STELLARTYPEPSTCE1, STELLARTYPEPSTCE2, 
                                    MASSPRECE1, MASSPRECE2, MASSPSTCE1, MASSPSTCE2, SEMIMAJORAXISPRECE, SEMIMAJORAXISPSTCE, CIRCULARIZATIONTIME],
                                names= ["CEs", "EVENTSCE","STELLARTYPEPRECE1", "STELLARTYPEPRECE2", "STELLARTYPEPSTCE1", "STELLARTYPEPSTCE2", 
                                    "MASSPRECE1", "MASSPRECE2", "MASSPSTCE1", "MASSPSTCE2", "SEMIMAJORAXISPRECE", "SEMIMAJORAXISPSTCE", "CIRCULARIZATIONTIME"])) 
hdu.append(SP_hdu)
hdu.append(MT_hdu)
hdu.append(CE_hdu)
hdu.close()

e = dt.datetime.now()
# Displays Time
current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
duration = e - s
print(" \n Finish time :", current_time)