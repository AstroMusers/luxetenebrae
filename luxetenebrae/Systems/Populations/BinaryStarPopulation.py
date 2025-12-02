from runpy import run_path
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
from luxetenebrae.Utils import utils as utils      # functions from utils.py
from luxetenebrae.Systems.BinarySystem import BinarySystem

compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')

class BinarySystemEvolution:
    """Class for detailed evolution of binary systems."""

    def __init__(self, mode, name=None):

        self.mode = mode
        self.name = name if name else f'BinarySystemEvolution_{mode}'
        self.binary_population = {}

        self.pathToData = os.path.join(compasRootDir, "luxetenebrae/runs", mode)
        self.pathToLogs = os.path.join(compasRootDir, "luxetenebrae/logs", mode)
        self.pathToFiles = os.path.join(compasRootDir, "luxetenebrae/files", mode)

    def read_data(self, raw_data=True, process_level='Purged', filetype='.h5'):
        """Read the binary population from the data files."""
        if raw_data:
            data_outputs = self._read_raw_data()
        else:
            data_outputs = self._read_processed_data(process_level, filetype)
        self.data_outputs = data_outputs
        return data_outputs
    
    def _read_raw_data(self):
        """Read raw data from detailed output files."""
        rawPath = os.path.join(self.pathToData, "Raw")
        runs= [x[0] for x in os.walk(rawPath) if 'Detailed_Output' in x[1]]  # Adjust the slice as needed
        if len(runs) == 0:
            raise ValueError(f"No raw data found in {rawPath}")
        print('runs:', runs)
        data_outputs= []
        for run in runs:  # Adjust the slice as needed
            out = [f for f in os.listdir(run + '/Detailed_Output') if ".h5" in f]
            for f in out:
                print("Reading file: ", run + '/Detailed_Output' + "/" + f)
                try:
                    data = h5.File(run + '/Detailed_Output' + "/" + f)
                    data_outputs.append(data)
                except:
                    continue
        self.data_outputs = data_outputs
        return data_outputs
    
    def _read_processed_data(self, process_level, filetype):
        """Read processed data from detailed output files."""
        processedPath = os.path.join(self.pathToData, f"Processed/{process_level}")
        runs= [x[0] for x in os.walk(processedPath) if 'Detailed_Output' in x[1]]  # Adjust the slice as needed
        if len(runs) == 0:
            raise ValueError(f"No processed data found in {processedPath}")
        print('runs:', runs)
        data_outputs= []
        for run in runs:  # Adjust the slice as needed
            out = [f for f in os.listdir(run + '/Detailed_Output') if f.endswith(filetype)]
            for f in out:
                print("Reading file: ", run + '/Detailed_Output' + "/" + f)
                try:
                    data = h5.File(run + '/Detailed_Output' + "/" + f)
                    data_outputs.append(data)
                except:
                    continue
        self.data_outputs = data_outputs
        return data_outputs
    
    def get_info(self):
        """Get information about the binary population."""
        info_dict = {
            'mode': self.mode,
            'name': self.name,
            'number_of_systems': len(self.data_outputs),
        }
        return info_dict

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

    def generate_binaries(self, requested_properties=None):
        """Generate a binary star population."""
        
        if self.data_outputs is None:
            raise ValueError("Data outputs not loaded. Please run read_data() first.")
        
        for Data in self.data_outputs:
            print("Processing system: ", Data.filename)
            seed = Data['SEED'][()]
            self.binary_population[seed] = BinarySystem(Data, requested_properties)

    