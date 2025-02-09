#### PARENT SCRIPT: secundus_processng_082824.py
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
import Luxetenebrae.scripts.calculations as calc      # functions from calculations.py
from astropy.io import fits
from astropy.table import Table

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
mode = ["Default/","Limited/", "Default_WD_Enabled/","Limited_WD_Enabled/" ]

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

    runs= [x[0] for x in os.walk(pathToData + '/Runs/' + mod) if "COMPAS" in x[0]]
    data_outputs = []
    for run in runs:
        out = [f for f in os.listdir(run) if ".h5" in f]
        try:
            data = h5.File(run + "/" + out[0])
            print("Reading file: ", run + "/" + out[0])
            data_outputs.append(data)
        except:
            continue

    # Keeps corresponding numerical values. 
    # SPs
    if not os.path.exists(pathToData + '/Files/' + mod +  str(s.strftime("%m.%d"))): 
        os.makedirs(pathToData + '/Files/' + mod  +  str(s.strftime("%m.%d")))
    directoryf = pathToData + '/Files/' + mod  +  str(s.strftime("%m.%d"))

    fit_filename = directoryf + "/secundus.fits"
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

    


    # f = open(pathToData + '/Files/' + mod  + str(s.strftime("%m.%d")) +  "_general_outputs.txt", "a")

    # f.writelines(["\n","\n Run :", start_time])
    print(f'Run :, {start_time}')

    i = 0

    for Data in data_outputs:

        c = dt.datetime.now()
        current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
        print("Batch (of 1 000 000 sys.)" + str(i) +  " start time :", current_time)

        SP = Data['BSE_System_Parameters']

        seedsSP = SP['SEED'][()]
        statusSP = SP['Evolution_Status'][()]

        stellarTypeZamsSP1   =  SP['Stellar_Type@ZAMS(1)'][()]
        stellarTypeZamsSP2   =  SP['Stellar_Type@ZAMS(2)'][()]
        stellarTypeSP1   =  SP['Stellar_Type(1)'][()]
        stellarTypeSP2   =  SP['Stellar_Type(2)'][()]  

        massZamsSP1 = SP['Mass@ZAMS(1)'][()] 
        massZamsSP2 = SP['Mass@ZAMS(2)'][()]

        semimajorAxisZamsSP = SP['SemiMajorAxis@ZAMS'][()] 

        SPs.extend(seedsSP)
        EVOLUTIONSTATSP.extend(statusSP)

        STELLARTYPEZAMSSP1.extend(stellarTypeZamsSP1)
        STELLARTYPEZAMSSP2.extend(stellarTypeZamsSP2)
        STELLARTYPESP1.extend(stellarTypeSP1)
        STELLARTYPESP2.extend(stellarTypeSP2)   

        MASSZAMSSP1.extend(massZamsSP1)
        MASSZAMSSP2.extend(massZamsSP2)

        SEMIMAJORAXISZAMSSP.extend(semimajorAxisZamsSP)

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

        semimajorAxispreMT = MT['SemiMajorAxis<MT'][()] 
        semimajorAxispstMT = MT['SemiMajorAxis>MT'][()]

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

        semimajorAxispreCE = CE['SemiMajorAxis<CE'][()] 
        semimajorAxispstCE = CE['SemiMajorAxis>CE'][()]

        timeCE = CE['Time'][()]
        
        radiuspreCE1 = CE['Radius(1)<CE'][()]
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



        # PreMT refers to the mass just before the first MT event, while pstMT refers to the mass just after the last MT event.
   

        c = dt.datetime.now()
        current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
        print("Batch (of 1 000 000 sys.)" + str(i) +  " end time :", current_time)
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
                                        MASSZAMSSP1,MASSZAMSSP2,SEMIMAJORAXISZAMSSP, EVOLUTIONSTATSP], 
                                    names=["SPs","STELLARTYPEZAMSSP1", "STELLARTYPEZAMSSP2", "STELLARTYPESP1", "STELLARTYPESP2",
                                        "MASSZAMSSP1", "MASSZAMSSP2",  "SEMIMAJORAXISZAMSSP", "EVOLUTIONSTATSP"]))
    
    MT_hdu = fits.BinTableHDU(Table(data=[MTs,STELLARTYPEPREMT1,STELLARTYPEPREMT2,STELLARTYPEPSTMT1,STELLARTYPEPSTMT2,
                                        MASSPREMT1,MASSPREMT2,MASSPSTMT1,MASSPSTMT2,SEMIMAJORAXISPREMT, SEMIMAJORAXISPSTMT,
                                        EVENTSMT,MASSTRANSFERHISTORY,RADIUSPREMT1, RADIUSPREMT2, RADIUSPSTMT1, RADIUSPSTMT2,
                                        TIMEPREMT, TIMEPSTMT, CEAFTERMT], 
                                    names=["MTs","STELLARTYPEPREMT1", "STELLARTYPEPREMT2", "STELLARTYPEPSTMT1", "STELLARTYPEPSTMT2",
                                        "MASSPREMT1", "MASSPREMT2", "MASSPSTMT1", "MASSPSTMT2", "SEMIMAJORAXISPREMT", "SEMIMAJORAXISPSTMT",
                                        "EVENTSMT","MASSTRANSFERHISTORY", "RADIUSPREMT1", "RADIUSPREMT2", "RADIUSPSTMT1", "RADIUSPSTMT2",
                                        "TIMEPREMT", "TIMEPSTMT","CEAFTERMT"]))
    CE_hdu = fits.BinTableHDU(Table(data=[CEs, EVENTSCE, STELLARTYPEPRECE1, STELLARTYPEPRECE2, STELLARTYPEPSTCE1, STELLARTYPEPSTCE2, 
                                        MASSPRECE1, MASSPRECE2, MASSPSTCE1, MASSPSTCE2, SEMIMAJORAXISPRECE, SEMIMAJORAXISPSTCE],
                                    names= ["CEs", "EVENTSCE","STELLARTYPEPRECE1", "STELLARTYPEPRECE2", "STELLARTYPEPSTCE1", "STELLARTYPEPSTCE2", 
                                        "MASSPRECE1", "MASSPRECE2", "MASSPSTCE1", "MASSPSTCE2", "SEMIMAJORAXISPRECE", "SEMIMAJORAXISPSTCE"])) 
    hdu.append(SP_hdu)
    hdu.append(MT_hdu)
    hdu.append(CE_hdu)
    hdu.close()

e = dt.datetime.now()
# Displays Time
current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
duration = e - s
print(" \n Finish time :", current_time)