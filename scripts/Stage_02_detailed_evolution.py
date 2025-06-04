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
    # matplotlib.rcParams['figure.figsize'] = (15,10)
    # matplotlib.rcParams['lines.markersize'] = 1
    # matplotlib.rcParams['font.size'] = 14
    # matplotlib.rcParams['legend.loc'] = "upper right"
    ps =  open(f'{pathToData}/Files/{mod}Detailed_Output/processed_systems.txt', 'a')
    psr = ps.readlines()  
    
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
    c = dt.datetime.now()
    current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
    if not os.path.exists(pathToData + '/Files/' + mod + '/Detailed_Output/' + str(s.strftime("%m.%d"))): 
        os.makedirs(pathToData + '/Files/' + mod + '/Detailed_Output/' +  str(s.strftime("%m.%d")))
    directoryf = pathToData + '/Files/' + mod  + '/Detailed_Output/' +   str(s.strftime("%m.%d"))


    


    # f = open(pathToData + '/Files/' + mod  + str(s.strftime("%m.%d")) +  "_general_outputs.txt", "a")

    # f.writelines(["\n","\n Run :", start_time])
    # print(f'Data outputs :, {data_outputs}')

    i = 0



    for Data in data_outputs[:1]:


        Data = h5.File(Data.filename, 'r')
        seed = Data['SEED'][:]

        if seed in psr:
            print(f"Seed {seed} already processed, skipping...")
            continue
        ps.write(f"{seed}\n")
        ps.flush()
        ps.close()

        print("Batch (of " + str(len(data_outputs)) + " sys.) " + str(i) +  " start time :", current_time)

        print(Data.keys())
        # print(Data['Stellar_Type(1)'][:])
        # print(Data['Stellar_Type(2)'][:])
        # print(len(Data['Stellar_Type(1)'][:]))
        stellar_type_1 = Data['Stellar_Type(1)'][:]
        stellar_type_2 = Data['Stellar_Type(2)'][:]
        radius_1 = Data['Radius(1)'][:]
        radius_2 = Data['Radius(2)'][:]
        time = Data['Time'][:]
        semimajoraxis = (Data['SemiMajorAxis'][:]*const.R_sun).to(u.au).value

  
        mask = np.array([is_ms_bh_pair(st1, st2) for st1, st2 in zip(stellar_type_1, stellar_type_2)])
        print(f"Number of MS + BH states in this system: {np.sum(mask)}")
        if np.sum(mask) != 0:

            fit_filename = directoryf + f"/secundus_{seed}.fits"
            hdu_pr = fits.PrimaryHDU()
            hdu_pr.writeto(fit_filename, overwrite=True)
            hdu = fits.open(fit_filename, mode='update')

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

            # Create a FITS table with the data
            col1 = fits.Column(name='Seed', format='K', array=seed)
            col2 = fits.Column(name='Time_BH_MS', format='D', unit='yr', array=time)
            col3 = fits.Column(name='Avg_Semimajor_Axis', format='D', unit='AU', array=semaj_average)
            col4 = fits.Column(name='Total_Lifetime', format='D', unit='yr', array=bhms_lifetime)
            # Create a binary table HDU
            hdu_1 = fits.BinTableHDU.from_columns([col1, col2, col3, col4])
            # Write the FITS file
            hdu.append(hdu_1)

            col1 = fits.Column(name='Time', format='D', unit='Myr', array=time)
            col2 = fits.Column(name='Stellar_Type_1', format='K', array=stellar_type_1)
            col3 = fits.Column(name='Stellar_Type_2', format='K', array=stellar_type_2)
            col4 = fits.Column(name='SemiMajorAxis', format='D', unit='AU', array=semimajoraxis)
            col5 = fits.Column(name='Radius1', format='D', unit='Rsun', array=radius_1)
            col6 = fits.Column(name='Radius2', format='D', unit='Rsun', array=radius_2)
            col7 = fits.Column(name='Time_Filtered', format='K', array=time_filtered)
            col8 = fits.Column(name='Stellar_Type_1_Filtered', format='K', array=stype1_filtered)
            col9 = fits.Column(name='Stellar_Type_2_Filtered', format='K', array=stype2_filtered)
            col10 = fits.Column(name='Radius1_Filtered', format='D', unit='Rsun', array=radius1_filtered)
            col11 = fits.Column(name='Radius2_Filtered', format='D', unit='Rsun', array=radius2_filtered)
            col12 = fits.Column(name='SemiMajorAxis_Filtered', format='D', unit='AU', array=semimajoraxis_filtered)

            
            # Create a binary table HDU
            hdu_2 = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12])
            # Write the FITS file
            hdu.append(hdu_2)
            hdu.flush()
            # Close the FITS file
            hdu.close()


            print(f"Data successfully written to {fit_filename}")

        i = i + 1




    #     SP = Data['BSE_System_Parameters']
    #     # print(SP.keys())



    #     seedsSP = SP['SEED'][()]
    #     statusSP = SP['Evolution_Status'][()]

    #     print(f"Batch (of {len(seedsSP)} sys.)" + str(i) +  " start time :", current_time)

    #     stellarTypeZamsSP1   =  SP['Stellar_Type@ZAMS(1)'][()]
    #     stellarTypeZamsSP2   =  SP['Stellar_Type@ZAMS(2)'][()]
    #     stellarTypeSP1   =  SP['Stellar_Type(1)'][()]
    #     stellarTypeSP2   =  SP['Stellar_Type(2)'][()]  

    #     massZamsSP1 = SP['Mass@ZAMS(1)'][()] 
    #     massZamsSP2 = SP['Mass@ZAMS(2)'][()]

    #     semimajorAxisZamsSP = SP['SemiMajorAxis@ZAMS'][()] # in AU
    #     eccentricityZamsSP = SP['Eccentricity@ZAMS'][()]

    #     radiusSP1 = massZamsSP1**(0.8) # in R_sun
    #     radiusSP2 = massZamsSP2**(0.8) # in R_sun
    #     print(radiusSP1)


        

    #     SPs.extend(seedsSP)
    #     EVOLUTIONSTATSP.extend(statusSP)

    #     STELLARTYPEZAMSSP1.extend(stellarTypeZamsSP1)
    #     STELLARTYPEZAMSSP2.extend(stellarTypeZamsSP2)
    #     STELLARTYPESP1.extend(stellarTypeSP1)
    #     STELLARTYPESP2.extend(stellarTypeSP2)   

    #     MASSZAMSSP1.extend(massZamsSP1)
    #     MASSZAMSSP2.extend(massZamsSP2)

    #     SEMIMAJORAXISZAMSSP.extend(semimajorAxisZamsSP)
    #     ECCENTRICITYZAMSSP.extend(eccentricityZamsSP)

    #     RADIUSSP1.extend(radiusSP1)
    #     RADIUSSP2.extend(radiusSP2)



    #     MT = Data['BSE_RLOF']

    #     seedsMT = MT['SEED'][()]
    #     eventsMT = MT['MT_Event_Counter'][()]

    #     stellarTypepreMT1   =  MT['Stellar_Type(1)<MT'][()]
    #     stellarTypepreMT2   =  MT['Stellar_Type(2)<MT'][()]  
    #     stellarTypepstMT1   =  MT['Stellar_Type(1)>MT'][()]
    #     stellarTypepstMT2   =  MT['Stellar_Type(2)>MT'][()] 

    #     masspreMT1 = MT['Mass(1)<MT'][()] 
    #     masspreMT2 = MT['Mass(2)<MT'][()] 
    #     masspstMT1 = MT['Mass(1)>MT'][()] 
    #     masspstMT2 = MT['Mass(2)>MT'][()]

    #     semimajorAxispreMT = MT['SemiMajorAxis<MT'][()] # in R_sun
    #     semimajorAxispstMT = MT['SemiMajorAxis>MT'][()]

    #     eccentricitypreMT = MT['Eccentricity<MT'][()]
    #     eccentricitypstMT = MT['Eccentricity>MT'][()]

    #     semimajorAxispreMT = (semimajorAxispreMT*const.R_sun).to(u.au).value
    #     semimajorAxispstMT = (semimajorAxispstMT*const.R_sun).to(u.au).value


    #     timepreMT = MT['Time<MT'][()]
    #     timepstMT = MT['Time>MT'][()]
        
    #     massTransferhistory = MT['CEE>MT'][()]

    #     radiuspreMT1 = MT['Radius(1)<MT'][()]
    #     radiuspstMT1 = MT['Radius(1)>MT'][()] 
    #     radiuspreMT2 = MT['Radius(2)<MT'][()]
    #     radiuspstMT2 = MT['Radius(2)>MT'][()]

        

    #     CEafterMT = MT['CEE>MT'][()]

    #     MTs.extend(seedsMT)
    #     EVENTSMT.extend(eventsMT)
      
    #     STELLARTYPEPREMT1.extend(stellarTypepreMT1)
    #     STELLARTYPEPREMT2.extend(stellarTypepreMT2) 
    #     STELLARTYPEPSTMT1.extend(stellarTypepstMT1)
    #     STELLARTYPEPSTMT2.extend(stellarTypepstMT2) 

    #     MASSPREMT1.extend(masspreMT1)
    #     MASSPREMT2.extend(masspreMT2)
    #     MASSPSTMT1.extend(masspstMT1)
    #     MASSPSTMT2.extend(masspstMT2)

    #     SEMIMAJORAXISPREMT.extend(semimajorAxispreMT)
    #     SEMIMAJORAXISPSTMT.extend(semimajorAxispstMT)

    #     ECCENTRICITYPREMT.extend(eccentricitypreMT)
    #     ECCENTRICITYPSTMT.extend(eccentricitypstMT)

    #     TIMEPREMT.extend(timepreMT)
    #     TIMEPSTMT.extend(timepstMT)

    #     MASSTRANSFERHISTORY.extend(massTransferhistory)
        

   
    #     CEAFTERMT.extend(CEafterMT)

    #     RADIUSPREMT1.extend(radiuspreMT1)
    #     RADIUSPREMT2.extend(radiuspreMT2)
    #     RADIUSPSTMT1.extend(radiuspstMT1)
    #     RADIUSPSTMT2.extend(radiuspstMT2)



    #     CE = Data['BSE_Common_Envelopes']

    #     seedsCE = CE['SEED'][()]
    #     eventsCE = CE['CE_Event_Counter'][()]
       
    #     stellarTypepreCE1   =  CE['Stellar_Type(1)<CE'][()]
    #     stellarTypepreCE2   =  CE['Stellar_Type(2)<CE'][()]  
    #     stellarTypepstCE1   =  CE['Stellar_Type(1)'][()]
    #     stellarTypepstCE2   =  CE['Stellar_Type(2)'][()] 
    #     # PreCE refers to the mass just before the first CE event, while pstCE refers to the mass just after the last CE event.
    #     masspreCE1 = CE['Mass(1)<CE'][()] 
    #     masspreCE2 = CE['Mass(2)<CE'][()] 
    #     masspstCE1 = CE['Mass(1)>CE'][()] 
    #     masspstCE2 = CE['Mass(2)>CE'][()]

    #     semimajorAxispreCE = CE['SemiMajorAxis<CE'][()] # in R_sun
    #     semimajorAxispstCE = CE['SemiMajorAxis>CE'][()]

    #     semimajorAxispreCE = (semimajorAxispreCE*const.R_sun).to(u.au).value
    #     semimajorAxispstCE = (semimajorAxispstCE*const.R_sun).to(u.au).value

    #     timeCE = CE['Time'][()]
    #     circulartizationTime = CE['Tau_Circ'][()]
        
    #     radiuspreCE1 = CE['Radius(1)<CE'][()] # in R_sun
    #     radiuspstCE1 = CE['Radius(1)>CE'][()] 
    #     radiuspreCE2 = CE['Radius(2)<CE'][()]
    #     radiuspstCE2 = CE['Radius(2)>CE'][()] 

    #     CEs.extend(seedsCE)
    #     EVENTSCE.extend(eventsCE)
   
    #     STELLARTYPEPRECE1.extend(stellarTypepreCE1)
    #     STELLARTYPEPRECE2.extend(stellarTypepreCE2) 
    #     STELLARTYPEPSTCE1.extend(stellarTypepstCE1)
    #     STELLARTYPEPSTCE2.extend(stellarTypepstCE2) 

    #     MASSPRECE1.extend(masspreCE1)
    #     MASSPRECE2.extend(masspreCE2)
    #     MASSPSTCE1.extend(masspstCE1)
    #     MASSPSTCE2.extend(masspstCE2)

    #     SEMIMAJORAXISPRECE.extend(semimajorAxispreCE)
    #     SEMIMAJORAXISPSTCE.extend(semimajorAxispstCE)

    #     CIRCULARIZATIONTIME.extend(circulartizationTime)


    #     # PreMT refers to the mass just before the first MT event, while pstMT refers to the mass just after the last MT event.
   

    #     c = dt.datetime.now()
    #     current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
    #     print("Batch " + str(i) +  " end time :", current_time)
    #     i = i+1
    # Data.close()

    # arrays_to_save = [ ]
    # filenames = []
    # checklist = [MTs,STELLARTYPEPREMT1,STELLARTYPEPREMT2,STELLARTYPEPSTMT1,STELLARTYPEPSTMT2,
    #                                     MASSPREMT1,MASSPREMT2,MASSPSTMT1,MASSPSTMT2,SEMIMAJORAXISPREMT, SEMIMAJORAXISPSTMT,
    #                                     EVENTSMT,MASSTRANSFERHISTORY,RADIUSPREMT1, RADIUSPREMT2, RADIUSPSTMT1, RADIUSPSTMT2,
    #                                     TIMEPREMT, TIMEPSTMT, CEAFTERMT]
    # for e in checklist:
    #     print(len(e))

    # SP_hdu = fits.BinTableHDU(Table(data=[SPs, STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,
    #                                     MASSZAMSSP1,MASSZAMSSP2,SEMIMAJORAXISZAMSSP, EVOLUTIONSTATSP, ECCENTRICITYZAMSSP, RADIUSSP1, RADIUSSP2], 
    #                                 names=["SPs","STELLARTYPEZAMSSP1", "STELLARTYPEZAMSSP2", "STELLARTYPESP1", "STELLARTYPESP2",
    #                                     "MASSZAMSSP1", "MASSZAMSSP2",  "SEMIMAJORAXISZAMSSP", "EVOLUTIONSTATSP", "ECCENTRICITYZAMSSP", "RADIUSSP1", "RADIUSSP2"]))
    
    # MT_hdu = fits.BinTableHDU(Table(data=[MTs,STELLARTYPEPREMT1,STELLARTYPEPREMT2,STELLARTYPEPSTMT1,STELLARTYPEPSTMT2,
    #                                     MASSPREMT1,MASSPREMT2,MASSPSTMT1,MASSPSTMT2,SEMIMAJORAXISPREMT, SEMIMAJORAXISPSTMT,
    #                                     EVENTSMT,MASSTRANSFERHISTORY,RADIUSPREMT1, RADIUSPREMT2, RADIUSPSTMT1, RADIUSPSTMT2,
    #                                     TIMEPREMT, TIMEPSTMT, CEAFTERMT, ECCENTRICITYPREMT, ECCENTRICITYPSTMT], 
    #                                 names=["MTs","STELLARTYPEPREMT1", "STELLARTYPEPREMT2", "STELLARTYPEPSTMT1", "STELLARTYPEPSTMT2",
    #                                     "MASSPREMT1", "MASSPREMT2", "MASSPSTMT1", "MASSPSTMT2", "SEMIMAJORAXISPREMT", "SEMIMAJORAXISPSTMT",
    #                                     "EVENTSMT","MASSTRANSFERHISTORY", "RADIUSPREMT1", "RADIUSPREMT2", "RADIUSPSTMT1", "RADIUSPSTMT2",
    #                                     "TIMEPREMT", "TIMEPSTMT","CEAFTERMT", "ECCENTRICITYPREMT", "ECCENTRICITYPSTMT"]))
    # CE_hdu = fits.BinTableHDU(Table(data=[CEs, EVENTSCE, STELLARTYPEPRECE1, STELLARTYPEPRECE2, STELLARTYPEPSTCE1, STELLARTYPEPSTCE2, 
    #                                     MASSPRECE1, MASSPRECE2, MASSPSTCE1, MASSPSTCE2, SEMIMAJORAXISPRECE, SEMIMAJORAXISPSTCE, CIRCULARIZATIONTIME],
    #                                 names= ["CEs", "EVENTSCE","STELLARTYPEPRECE1", "STELLARTYPEPRECE2", "STELLARTYPEPSTCE1", "STELLARTYPEPSTCE2", 
    #                                     "MASSPRECE1", "MASSPRECE2", "MASSPSTCE1", "MASSPSTCE2", "SEMIMAJORAXISPRECE", "SEMIMAJORAXISPSTCE", "CIRCULARIZATIONTIME"])) 
    # hdu.append(SP_hdu)
    # hdu.append(MT_hdu)
    # hdu.append(CE_hdu)
    # hdu.close()

e = dt.datetime.now()
# Displays Time
current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
duration = e - s
print(" \n Finish time :", current_time)