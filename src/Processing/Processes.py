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
from src.Utils import utils as utils      # functions from utils.py
import astropy.units as u
from astropy import constants as const
# Import COMPAS specific scripts
import astropy.io.fits as fits
from astropy.table import Table
from src.Systems.Populations.BinaryStarPopulation import BinarySystemEvolution
import src.Utils.calculations as calc

class Process:
    """Class to process the COMPAS output data, including purging non-BH systems and extracting system data."""
    
    def __init__(self, binary_system_evolution: BinarySystemEvolution):

        self.bse = binary_system_evolution
        self.name = f"Process_{self.bse.mode}"

    def purge_nonbh_systems(self, delete_after_purge=True, with_logging=True):
        """Purge non-black hole systems from the detailed output files."""
        # Displays Time
        # Choose the mode to process
        now = dt.datetime.now()

        # Create output directory if it doesn't exist
        os.makedirs(self.bse.pathToFiles, exist_ok=True)

        logger = logging.getLogger(self.name)
        logger.setLevel(logging.INFO)

        fh = None
        if with_logging:
            os.makedirs(self.bse.pathToLogs, exist_ok=True)
            log_filename = os.path.join(self.bse.pathToLogs, f"purge_{now.strftime('%Y%m%d_%H%M%S')}.log")
            fh = logging.FileHandler(log_filename, mode='w')
            fh.setLevel(logging.INFO)
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            fh.setFormatter(formatter)
            # Prevent duplicate handlers if method called multiple times
            if not any(isinstance(h, logging.FileHandler) and h.baseFilename == fh.baseFilename for h in logger.handlers):
                logger.addHandler(fh)
            logger.propagate = False  # keep logs out of root logger

            # Optional: redirect prints to logger (comment out in Jupyter if undesired)
            sys.stdout = utils.StreamToLogger(logger, logging.INFO)
            sys.stderr = utils.StreamToLogger(logger, logging.ERROR)
            logger.info('Script purge started')

        # Import COMPAS specific scripts

        # print("COMPAS_ROOT_DIR:", compasRootDir)
        # print(sys.path)



        s = dt.datetime.now()
        start_time = s.strftime("%d%m%y") + "_" + s.strftime('%H%M')
        print("Start time :", start_time)
        n_systems = len(self.bse.data_outputs)
        print(f"Total number of systems to process: {n_systems}")
        purgePath = os.path.join(self.bse.pathToData, f"Processed/Purged")
        os.makedirs(purgePath, exist_ok=True)


        i = j = k = l = b = s = er = 0
        # Loop through the data outputs and check for black holes and MS + BH pairs. If not found, delete the detailed output file. 
        print(f"Total number of systems to process: {len(self.bse.detailed_data_outputs)}")
        print("Starting the purge...")

        for Data in self.bse.detailed_data_outputs:

            try:

                k += 1
                Data = h5.File(Data.filename, 'r')
                seed = Data['SEED'][()]

                print(f"System of {k} (of {len(self.bse.detailed_data_outputs)} sys.)")



                stellar_type_1 = Data['Stellar_Type(1)'][:]
                stellar_type_2 = Data['Stellar_Type(2)'][:]
                BH1 = [x for x in stellar_type_1 if x == 14]
                BH2 = [x for x in stellar_type_2 if x == 14]

                print(f"System of {k} (of {len(self.bse.detailed_data_outputs)} sys.) - Number of primary BH: {len(BH1)}, Number of secondary BH: {len(BH2)}")

                # Check if there are any BHs in the system
                if (len(BH1) + len(BH2)) == 0:
                    if delete_after_purge:
                        print(f"No BHs in this system, deleting {Data.filename}")
                        os.remove(f'{Data.filename}')
                    else:
                        print(f"No BHs in this system")
                    i += 1
                    continue

                # Check if there are any MS + BH pairs in the system
                mask = np.array([utils.is_ms_bh_pair(st1, st2) for st1, st2 in zip(stellar_type_1, stellar_type_2)])
                print(f"Number of MS + BH states in this system: {np.sum(mask)}")
                if np.sum(mask) == 0:
                    if delete_after_purge:
                        print(f"No BH-MS states occur in this batch, deleting {Data.filename}")
                        os.remove(f'{Data.filename}')
                    else:
                        print(f"No BH-MS states occur in this batch")
                    j += 1
                    continue

                # Check if there are any Bh + BH pairs in the system
                mask_bh = np.array([utils.is_bh_bh_pair(st1, st2) for st1, st2 in zip(stellar_type_1, stellar_type_2)])
                print(f"Number of BH + BH states in this system: {np.sum(mask_bh)}")
                if np.sum(mask_bh) > 0:
                    print(f"System {seed[0]} has BH + BH pairs, keeping it.")
                    b = b + 1
                
                # Save the seed of the system
                print(f"System {seed[0]} survived the purge.")
                # Move the file to the purged directory
                unique_filename = f"{seed[0]}_{os.path.basename(Data.filename)}"
                dest_file = os.path.join(purgePath, unique_filename)
                print(f"Moving file to {dest_file}")
                os.rename(Data.filename, dest_file)
                l = l + 1

            except Exception as e:
                print(f"Error processing system {k}: {e}")
                logging.error(f"Error processing system {k}: {e}")
                er = er + 1
                continue
                    

        print(f"Total number of systems read: {k}")
        print(f"Total number of systems skipped: {s}")
        print(f"Total number of systems without BHs: {i}")
        print(f"Total number of systems without BH-MS pairs: {j}")           
        print(f"Total number of systems that survived the purge: {l}")
        print(f"Total number of BH + BH pairs in the systems: {b}")
        print(f"Total number of errors during processing: {er}")

        if with_logging:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            logging.shutdown()

        results = [f"Total number of systems: {k} Total number of systems without BHs: {i} Total number of systems without BH-MS pairs: {j} Total number of systems that survived the purge: {l} Total number of BH + BH pairs in the systems: {b} survivors file created at: {self.bse.pathToFiles}/purged_{start_time}.txt"
                ]

        return results

    def extract_detailed_system_data(self, process_from='Purged', with_logging=True):
        
        start_time = dt.datetime.now()

        # Create output directory if it doesn't exist
        os.makedirs(self.bse.pathToFiles, exist_ok=True)

        logger = logging.getLogger(self.name)
        logger.setLevel(logging.INFO)

        fh = None
        if with_logging:
            os.makedirs(self.bse.pathToLogs, exist_ok=True)
            log_filename = os.path.join(self.bse.pathToLogs, f"extract_detailed_system_data_{start_time.strftime('%Y%m%d_%H%M%S')}.log")
            fh = logging.FileHandler(log_filename, mode='w')
            fh.setLevel(logging.INFO)
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            fh.setFormatter(formatter)
            # Prevent duplicate handlers if method called multiple times
            if not any(isinstance(h, logging.FileHandler) and h.baseFilename == fh.baseFilename for h in logger.handlers):
                logger.addHandler(fh)
            logger.propagate = False  # keep logs out of root logger

            # Optional: redirect prints to logger (comment out in Jupyter if undesired)
            sys.stdout = utils.StreamToLogger(logger, logging.INFO)
            sys.stderr = utils.StreamToLogger(logger, logging.ERROR)
            logger.info('data extraction started')

        # Check required output data exists
        # if 'Purged' in process_from:
        detailed_data = self.bse.get_processed_data(level=process_from)
        print(f"Extracting system data from {len(detailed_data)} purged systems.")
        # else:
        #     detailed_data = self.bse.detailed_data_outputs
        #     print(f"Extracting system data from raw {len(detailed_data)} systems.")

        base_data = self.bse.data_outputs
        
        start_time_str = start_time.strftime("%d%m%y") + "_" + start_time.strftime('%H%M')
        print("Start time :", start_time_str)
        n_systems = len(detailed_data)
        print(f"Total number of systems to process: {n_systems}")

        i = 0
        j = 0
        k = 0
        l = 0
        SP_dict = {}
        MT_dict = {}
        CE_dict = {}

        survivor_seeds = []
        survivors = []
        for Data in detailed_data:
            Data = h5.File(Data.filename, 'r')
            seeds = Data['SEED'][()]
            seed = seeds[0]  # Assuming we want to process the first seed in the file
            survivor_seeds.append(seed)
            survivors.append(Data)

        print(f"Total number of systems that survived the purge: {len(survivor_seeds)}")

        for Data in base_data:

            Data = h5.File(Data.filename, 'r')
            SP_d, MT_d, CE_d = utils.extract_masked_systems_by_seed(Data, survivor_seeds) 
            SP_dict.update(SP_d)
            MT_dict.update(MT_d)
            CE_dict.update(CE_d)

        for Data, seed in zip(survivors, survivor_seeds):


            fit_filename = os.path.join(self.bse.pathToFiles, f"{process_from}_processed__detailed_outputs_{seed}.fits") 
            hdu_pr = fits.PrimaryHDU()
            hdu_pr.writeto(fit_filename, overwrite=True)
            hdu = fits.open(fit_filename, mode='update')

            seeds = Data['SEED'][()]
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

            table_one = fits.BinTableHDU(Table(data=[seeds, time, stellar_type_1, stellar_type_2, semimajoraxis, radius_1, radius_2,
                                                roche_lobe_1, roche_lobe_2, mass_1, mass_2, mass_zams_1, mass_zams_2, mt_history, eccentricity],
                                            names=['Seed', 'Time', 'Stellar_Type_1', 'Stellar_Type_2', 'SemiMajorAxis',
                                                    'Radius(1)', 'Radius(2)', 'RocheLobe(1)', 'RocheLobe(2)', 'Mass(1)', 'Mass(2)',
                                                    'Mass@ZAMS(1)', 'Mass@ZAMS(2)', 'MT_History', 'Eccentricity'],
                                            units=[' ', 'Myr', ' ', ' ', 'AU', 'Rsun', 'Rsun', 'Rsun', 'Rsun', 'Msun', 'Msun', 'Msun', 'Msun', ' ', ' ']))
            hdu.append(table_one)
            hdu.flush()

            seedSP = [SP_dict[seed]['SEED'][()]]  # Convert to list for consistency
            statusSP = [SP_dict[seed]['Evolution_Status'][()]]
            stellarTypeZamsSP1 = [SP_dict[seed]['Stellar_Type@ZAMS(1)'][()]]
            stellarTypeZamsSP2 = [SP_dict[seed]['Stellar_Type@ZAMS(2)'][()]]
            stellarTypeSP1 = [SP_dict[seed]['Stellar_Type(1)'][()]]
            stellarTypeSP2 = [SP_dict[seed]['Stellar_Type(2)'][()]]
            massZamsSP1 = [SP_dict[seed]['Mass@ZAMS(1)'][()]]
            massZamsSP2 = [SP_dict[seed]['Mass@ZAMS(2)'][()]]
            semimajorAxisZamsSP = [SP_dict[seed]['SemiMajorAxis@ZAMS'][()]]  # in AU
            eccentricityZamsSP = [SP_dict[seed]['Eccentricity@ZAMS'][()]]
            mergersSP = [SP_dict[seed]['Merger'][()]]

            table_two = fits.BinTableHDU(Table(data=[seedSP, statusSP, stellarTypeZamsSP1, stellarTypeZamsSP2, stellarTypeSP1, stellarTypeSP2,
                                                massZamsSP1, massZamsSP2, semimajorAxisZamsSP, eccentricityZamsSP, mergersSP],
                                            names=['SEED', 'Evolution_Status', 'Stellar_Type@ZAMS(1)', 'Stellar_Type@ZAMS(2)',
                                                    'Stellar_Type(1)', 'Stellar_Type(2)', 'Mass@ZAMS(1)', 'Mass@ZAMS(2)',
                                                    'SemiMajorAxis@ZAMS', 'Eccentricity@ZAMS', 'Merger'],
                                            units=[' ', ' ', ' ', ' ', ' ', ' ', 'Msun', 'Msun', 'AU', '', '']))
            hdu.append(table_two)
            hdu.flush()

            if seed not in MT_dict:
                print(f"Seed {seed} not in MT_dict, third table not created.")
                j = j + 1
            
            else:

                seedMT = [MT_dict[seed]['SEED_MT'][()]]
                eventsMT = [MT_dict[seed]['MT_Event_Counter'][()]]
                stellarTypepreMT1 = [MT_dict[seed]['Stellar_Type(1)<MT'][()]]
                stellarTypepreMT2 = [MT_dict[seed]['Stellar_Type(2)<MT'][()]]
                stellarTypepstMT1 = [MT_dict[seed]['Stellar_Type(1)>MT'][()]]
                stellarTypepstMT2 = [MT_dict[seed]['Stellar_Type(2)>MT'][()]]
                masspreMT1 = [MT_dict[seed]['Mass(1)<MT'][()]]
                masspreMT2 = [MT_dict[seed]['Mass(2)<MT'][()]]
                masspstMT1 = [MT_dict[seed]['Mass(1)>MT'][()]]
                masspstMT2 = [MT_dict[seed]['Mass(2)>MT'][()]]
                semimajorAxispreMT = [MT_dict[seed]['SemiMajorAxis<MT'][()]]  # in AU
                semimajorAxispstMT = [MT_dict[seed]['SemiMajorAxis>MT'][()]]
                eccentricitypreMT = [MT_dict[seed]['Eccentricity<MT'][()]]
                eccentricitypstMT = [MT_dict[seed]['Eccentricity>MT'][()]]

                timepreMT = [MT_dict[seed]['Time<MT'][()]]
                timepstMT = [MT_dict[seed]['Time>MT'][()]]
                radiuspreMT1 = [MT_dict[seed]['Radius(1)<MT'][()]]
                radiuspreMT2 = [MT_dict[seed]['Radius(2)<MT'][()]]
                radiuspstMT1 = [MT_dict[seed]['Radius(1)>MT'][()]]
                radiuspstMT2 = [MT_dict[seed]['Radius(2)>MT'][()]]
                CEafterMT = [MT_dict[seed]['CE_History'][()]]

                table_three = fits.BinTableHDU(Table(data=[seedMT, eventsMT, stellarTypepreMT1, stellarTypepreMT2, stellarTypepstMT1, stellarTypepstMT2,
                                                    masspreMT1, masspreMT2, masspstMT1, masspstMT2, semimajorAxispreMT, semimajorAxispstMT,
                                                    eccentricitypreMT, eccentricitypstMT, timepreMT, timepstMT,
                                                    radiuspreMT1, radiuspreMT2, radiuspstMT1, radiuspstMT2, CEafterMT],
                                                names=['SEED', 'MT_Event_Counter', 'Stellar_Type(1)<MT', 'Stellar_Type(2)<MT',
                                                        'Stellar_Type(1)>MT', 'Stellar_Type(2)>MT', 'Mass(1)<MT', 'Mass(2)<MT',
                                                        'Mass(1)>MT', 'Mass(2)>MT', 'SemiMajorAxis<MT', 'SemiMajorAxis>MT',
                                                        'Eccentricity<MT', 'Eccentricity>MT', 'Time<MT', 'Time>MT',
                                                        'Radius(1)<MT', 'Radius(2)<MT', 'Radius(1)>MT', 'Radius(2)>MT',
                                                        'CE_History'],
                                                units=[' ', '', '', '', '', '', 'Msun', 'Msun', 'Msun', 'Msun', 'AU', 'AU',
                                                        '', '', 'Myr', 'Myr', 'Rsun', 'Rsun', 'Rsun', 'Rsun', '']))
                hdu.append(table_three)
                hdu.flush()  

            if seed not in CE_dict:
                print(f"Seed {seed} not in CE_dict, fourth table not created.")
                k = k + 1
            else:

                seedCE = [CE_dict[seed]['SEED_CE'][()]]
                eventsCE = [CE_dict[seed]['CE_Event_Counter'][()]]
                stellarTypepreCE1 = [CE_dict[seed]['Stellar_Type(1)<CE'][()]]
                stellarTypepreCE2 = [CE_dict[seed]['Stellar_Type(2)<CE'][()]]
                stellarTypepstCE1 = [CE_dict[seed]['Stellar_Type(1)'][()]]
                stellarTypepstCE2 = [CE_dict[seed]['Stellar_Type(2)'][()]]
                masspreCE1 = [CE_dict[seed]['Mass(1)<CE'][()]]
                masspreCE2 = [CE_dict[seed]['Mass(2)<CE'][()]]
                masspstCE1 = [CE_dict[seed]['Mass(1)>CE'][()]]
                masspstCE2 = [CE_dict[seed]['Mass(2)>CE'][()]]
                semimajorAxispreCE = [CE_dict[seed]['SemiMajorAxis<CE'][()]]  # in au
                semimajorAxispstCE = [CE_dict[seed]['SemiMajorAxis>CE'][()]]
                timeCE = [CE_dict[seed]['Time_CE'][()]]
                circularizationTime = [CE_dict[seed]['Tau_Circ'][()]]
                radiuspreCE1 = [CE_dict[seed]['Radius(1)<CE'][()]]
                radiuspreCE2 = [CE_dict[seed]['Radius(2)<CE'][()]]
                radiuspstCE1 = [CE_dict[seed]['Radius(1)>CE'][()]]
                radiuspstCE2 = [CE_dict[seed]['Radius(2)>CE'][()]]

                table_four = fits.BinTableHDU(Table(data=[seedCE, eventsCE, stellarTypepreCE1, stellarTypepreCE2, stellarTypepstCE1, stellarTypepstCE2,
                                                    masspreCE1, masspreCE2, masspstCE1, masspstCE2, semimajorAxispreCE, semimajorAxispstCE,
                                                    timeCE, circularizationTime,
                                                    radiuspreCE1, radiuspreCE2, radiuspstCE1, radiuspstCE2],
                                                names=['SEED', 'CE_Event_Counter', 'Stellar_Type(1)<CE', 'Stellar_Type(2)<CE',
                                                        'Stellar_Type(1)>CE', 'Stellar_Type(2)>CE', 'Mass(1)<CE', 'Mass(2)<CE',
                                                        'Mass(1)>CE', 'Mass(2)>CE', 'SemiMajorAxis<CE', 'SemiMajorAxis>CE',
                                                        'Time_CE', 'Tau_Circ',
                                                        'Radius(1)<CE', 'Radius(2)<CE', 'Radius(1)>CE', 'Radius(2)>CE'],
                                                units=[' ', '', '', '', '', '', 'Msun', 'Msun', 'Msun', 'Msun', 'AU', 'AU',
                                                        'Myr', 'Myr', 'Rsun', 'Rsun', 'Rsun', 'Rsun']))
                hdu.append(table_four)
                hdu.flush()

            print(f"System of {l} (of {len(survivor_seeds)} sys.)")
            # Close the FITS file
            hdu.close()


            print(f"Data successfully written to {fit_filename}")

            l = l + 1

        print(f"Total number of systems with BHs: {l}")
        print(f"Total number of systems without MT data: {j}")
        print(f"Total number of systems without CE data: {k}")
        print(f"Total number of systems without CE data: {k}")
        end_time = dt.datetime.now()
        # Displays Time
        end_time_str = end_time.strftime("%d%m%y") + "_" + end_time.strftime('%H%M')
        duration = end_time.timestamp() - start_time.timestamp()
        print(" \n Finish time :", end_time_str, f"Duration: {duration}")

        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        logging.shutdown()

        results = [f"Total number of systems with BHs: {l}",
                f"Total number of systems without MT data: {j}",
                f"Total number of systems without CE data: {k}",
                f"Output files created in: {self.bse.pathToFiles}"]
        
        return results
    
    def extract_base_system_data(self, with_logging=True):
        """Extracts system data from COMPAS detailed output files and saves it to a FITS file."""
        start_time = dt.datetime.now()

        # Create output directory if it doesn't exist
        os.makedirs(self.bse.pathToFiles, exist_ok=True)

        logger = logging.getLogger(self.name)
        logger.setLevel(logging.INFO)

        fh = None
        if with_logging:
            os.makedirs(self.bse.pathToLogs, exist_ok=True)
            log_filename = os.path.join(self.bse.pathToLogs, f"extract_base_system_data_{start_time.strftime('%Y%m%d_%H%M%S')}.log")
            fh = logging.FileHandler(log_filename, mode='w')
            fh.setLevel(logging.INFO)
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            fh.setFormatter(formatter)
            # Prevent duplicate handlers if method called multiple times
            if not any(isinstance(h, logging.FileHandler) and h.baseFilename == fh.baseFilename for h in logger.handlers):
                logger.addHandler(fh)
            logger.propagate = False  # keep logs out of root logger

            # Optional: redirect prints to logger (comment out in Jupyter if undesired)
            sys.stdout = utils.StreamToLogger(logger, logging.INFO)
            sys.stderr = utils.StreamToLogger(logger, logging.ERROR)
            logger.info('data extraction started')

        base_data = self.bse.data_outputs


        fit_filename = os.path.join(self.bse.pathToFiles, f"processed_base_outputs.fits")
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

        print(f'Data outputs :, {base_data}')

        i = 0

        for Data in base_data:

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

        end_time = dt.datetime.now()
        # Displays Time
        end_time_str = end_time.strftime("%d%m%y") + "_" + end_time.strftime('%H%M')
        duration = end_time.timestamp() - start_time.timestamp()
        print(" \n Finish time :", end_time_str, f"Duration: {duration}")

        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        logging.shutdown()

    def create_base_data_masks(self, with_logging=True):
        """Create data masks for heatmaps."""

        start_time = dt.datetime.now()

        # Create output directory if it doesn't exist
        os.makedirs(self.bse.pathToFiles, exist_ok=True)

        logger = logging.getLogger(self.name)
        logger.setLevel(logging.INFO)

        fh = None
        if with_logging:
            os.makedirs(self.bse.pathToLogs, exist_ok=True)
            log_filename = os.path.join(self.bse.pathToLogs, f"mask_base_data_{start_time.strftime('%Y%m%d_%H%M%S')}.log")
            fh = logging.FileHandler(log_filename, mode='w')
            fh.setLevel(logging.INFO)
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            fh.setFormatter(formatter)
            # Prevent duplicate handlers if method called multiple times
            if not any(isinstance(h, logging.FileHandler) and h.baseFilename == fh.baseFilename for h in logger.handlers):
                logger.addHandler(fh)
            logger.propagate = False  # keep logs out of root logger

            # Optional: redirect prints to logger (comment out in Jupyter if undesired)
            sys.stdout = utils.StreamToLogger(logger, logging.INFO)
            sys.stderr = utils.StreamToLogger(logger, logging.ERROR)
            logger.info('data extraction started')


        start_time_str = start_time.strftime("%d%m%y") + "_" + start_time.strftime('%H%M')
        print("Start time :", start_time_str)
        # Choose an output hdf5 file to work with

        fit_filename = os.path.join(self.bse.pathToFiles, f"masks_base_outputs.fits")
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

        fits_data = os.path.join(self.bse.pathToFiles, f"processed_base_outputs.fits")

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

        MASKMTinSP =  np.isin(SP['SPs'], MT['MTs'])

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

        MASKPREMTnegCE =  np.isin(MT['MTs']*MASKFIRSTMT, CE['CEs'], invert=True)
        MASKPREMTposCE =  np.isin(MT['MTs']*MASKFIRSTMT, CE['CEs'])
        MASKPSTMTnegCE =  np.isin(MT['MTs']*MASKLASTMT, CE['CEs'], invert=True)
        MASKPSTMTposCE =  np.isin(MT['MTs']*MASKLASTMT, CE['CEs'])
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


        end_time = dt.datetime.now()
        # Displays Time
        end_time_str = end_time.strftime("%d%m%y") + "_" + end_time.strftime('%H%M')
        duration = end_time.timestamp() - start_time.timestamp()
        print(" \n Finish time :", end_time_str, f"Duration: {duration}")

        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        logging.shutdown()