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

class Process(BinarySystemEvolution):
    """
    A class to represent the processing of binary star systems in the luxetenebrae package.
    Inherits from BinarySystemEvolution.
    Public attributes: None (inherits from BinarySystemEvolution)
    Internal (prefixed with _): None (inherits from BinarySystemEvolution)
    """
    def __init__(self, mode, name=None):
        super().__init__(mode=mode, name=name)
        self._set_directory_paths(mode)
        self.pathToData, self.pathToLogs, self.pathToFiles = self._set_directory_paths(mode)

    def _set_directory_paths(mode):
        """Set directory paths based on the mode."""
        compasRootDir = os.environ['COMPAS_ROOT_DIR']
        pathToData = os.path.join(compasRootDir, "../luxetenebrae/runs", mode)
        pathToLogs = os.path.join(compasRootDir, "../luxetenebrae/logs", mode)
        pathToFiles = os.path.join(compasRootDir, "../luxetenebrae/files", mode)
        return pathToData, pathToLogs, pathToFiles

    def purge_nonbh_systems(self, delete_after_purge=True, with_logging=True):
        """Purge non-black hole systems from the detailed output files."""
        # Displays Time
        # Choose the mode to process
        now = dt.datetime.now()

        # Create output directory if it doesn't exist
        os.makedirs(self.pathToFiles, exist_ok=True)

        logger = logging.getLogger(self.name)
        logger.setLevel(logging.INFO)

        fh = None
        if with_logging:
            os.makedirs(self.pathToLogs, exist_ok=True)
            log_filename = os.path.join(self.pathToLogs, f"purge_{now.strftime('%Y%m%d_%H%M%S')}.log")
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
        n_systems = len(self.data_outputs)
        print(f"Total number of systems to process: {n_systems}")
        purgePath = os.path.join(self.pathToData, f"Processed/Purged")
        os.makedirs(purgePath, exist_ok=True)


        i = j = k = l = b = s = er = 0
        # Loop through the data outputs and check for black holes and MS + BH pairs. If not found, delete the detailed output file. 
        print(f"Total number of systems to process: {len(self.data_outputs)}")
        print("Starting the purge...")

        for Data in self.data_outputs:

            try:

                k += 1
                Data = h5.File(Data.filename, 'r')
                seed = Data['SEED'][()]

                print(f"System of {k} (of {len(self.data_outputs)} sys.)")



                stellar_type_1 = Data['Stellar_Type(1)'][:]
                stellar_type_2 = Data['Stellar_Type(2)'][:]
                BH1 = [x for x in stellar_type_1 if x == 14]
                BH2 = [x for x in stellar_type_2 if x == 14]

                print(f"System of {k} (of {len(self.data_outputs)} sys.) - Number of primary BH: {len(BH1)}, Number of secondary BH: {len(BH2)}")

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
                dest_file = os.path.join(purgePath, os.path.basename(Data.filename))
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

        results = [f"Total number of systems: {k} Total number of systems without BHs: {i} Total number of systems without BH-MS pairs: {j} Total number of systems that survived the purge: {l} Total number of BH + BH pairs in the systems: {b} survivors file created at: {self.pathToFiles}/purged_{start_time}.txt"
                ]

        return results

    def extract_system_data(self, proccess_from='Purged', with_logging=True):
        
        now = dt.datetime.now()

        # Create output directory if it doesn't exist
        os.makedirs(self.pathToFiles, exist_ok=True)

        logger = logging.getLogger(self.name)
        logger.setLevel(logging.INFO)

        fh = None
        if with_logging:
            os.makedirs(self.pathToLogs, exist_ok=True)
            log_filename = os.path.join(self.pathToLogs, f"extract_system_data_{now.strftime('%Y%m%d_%H%M%S')}.log")
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

        # Import COMPAS specific scripts

        # print("COMPAS_ROOT_DIR:", compasRootDir)
        # print(sys.path)



        s = dt.datetime.now()
        start_time = s.strftime("%d%m%y") + "_" + s.strftime('%H%M')
        print("Start time :", start_time)
        n_systems = len(self.data_outputs)
        print(f"Total number of systems to process: {n_systems}")

        if proccess_from == 'Purged':
            
        purgePath = os.path.join(self.pathToData, f"Processed/Purged")
        os.makedirs(purgePath, exist_ok=True)


        i = j = k = l = b = s = er = 0
        # Loop through the data outputs and check for black holes and MS + BH pairs. If not found, delete the detailed output file. 
        print(f"Total number of systems to process: {len(self.data_outputs)}")
        print("Starting the purge...")


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

        ps = open(f'{compasRootDir}/luxetenebrae/files/survivors.txt', 'r')
        survivors = ps.read().splitlines()
        runs= [x[0] for x in os.walk(pathToData) if "Detailed_Output" in x[1]]
        print('runs:', runs)
        data_outputs_1 = []
        data_outputs_2 = []
        for run in runs:  # Adjust the slice as needed
            out1 = [f for f in os.listdir(run + '/Detailed_Output') if ".h5" in f]
            out2 = [f for f in os.listdir(run) if ".h5" in f]
            for f in out1:
                print("Reading file: ", run + '/Detailed_Output' + "/" + f)
                try:
                    data = h5.File(run + '/Detailed_Output' + "/" + f)
                    data_outputs_1.append(data)
                except:
                    continue
            for f in out2:
                print("Reading file: ", run  + "/" + f)
                data = h5.File( run + "/" + f, 'r')
                data_outputs_2.append(data)
        # Keeps corresponding numerical values. 
        # SPs
        c = dt.datetime.now()
        current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')

        i = 0
        j = 0
        k = 0
        l = 0
        SP_dict = {}
        MT_dict = {}
        CE_dict = {}

        for Data in data_outputs_2:

            Data = h5.File(Data.filename, 'r')

            SP_d, MT_d, CE_d = utils.extract_masked_systems_by_seed(Data, survivors) 
            SP_dict.update(SP_d)
            MT_dict.update(MT_d)
            CE_dict.update(CE_d)

        for Data in data_outputs_1:

            Data = h5.File(Data.filename, 'r')
            seeds = Data['SEED'][()]
            seed = seeds[0]  # Assuming we want to process the first seed in the file

            # Check if the seed exists in the survivors list
            if str(seed) not in survivors:
                print(f"Seed {seed} not in survivors, skipping...")
                i = i + 1
                continue

            fit_filename = pathToFiles + f"/stage_02DE_outputs_{seed}.fits" 
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

            print(f"System of {l} (of {len(data_outputs_1)} sys.)")
            # Close the FITS file
            hdu.close()


            print(f"Data successfully written to {fit_filename}")

            l = l + 1
        print(f"Total number of systems with BHs: {l}")
        print(f"Total number of systems without MT data: {j}")
        print(f"Total number of systems without CE data: {k}")
        print(f"Total number of systems without CE data: {k}")
        e = dt.datetime.now()
        # Displays Time
        current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
        duration = e - s
        print(" \n Finish time :", current_time)

        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        logging.shutdown()

        results = [f"Total number of systems with BHs: {l}",
                f"Total number of systems without MT data: {j}",
                f"Total number of systems without CE data: {k}",
                f"Output files created in: {pathToFiles}"]
        return results
