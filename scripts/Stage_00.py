#PARENT SCRIPT: my_runSubmit_Default_a
import numpy as np
import sys
import os
from subprocess import call
import yaml
import argparse
import warnings
import logging
import datetime as dt


# The config file (my_compasConfigDefault.yaml) is used to set the options for COMPAS. You can change the options in the config file to customize the runs. 
# In config file, number of systems each run is set to 1000. Therefore, total number of systems generated with this script will be 1000 * number_of_runs.
# If you want to change the number of runs, change the number_of_runs variable.
# If you want to change the number of systems per run, change the value in the config file  under numericalChoices -> '--number-of-systems'.
# Detailed outputs is enabled by default in the config file, so you will get detailed outputs for each run. Change the option from 'True' to 'False' in the config file if you want to disable detailed outputs.
# The output directory is set to the current working directory by default. You can change the output directory by setting the environment variable COMPAS_LOGS_OUTPUT_DIR_PATH or by changing the value in the config file under stringChoices -> '--output-path'.
# Depending on your choices in the config file, you can keep track of important changes or categorise the runs by changing the mode variable below. In this way runs with different configurations can be stored in different directories.
number_of_runs = 5
mode = "WD_Enabled_Detailed"  # Change this to the desired 'mode'. We use 'WD_Enabled_Detailed' because we want to enable detailed outputs and we want to include white dwarfs in the runs.
for i in range(number_of_runs):
    # Check if we are using python 3
    python_version = sys.version_info[0]
    print("python_version =", python_version)

    HERE = os.path.dirname(__file__)
    DEFAULT_CONFIG_FILE = os.path.join(os.getcwd(), 'luxetenebrae/config/my_compasConfigDefault.yaml')
    REPO_ROOT = os.path.abspath(os.path.join(HERE, "../../"))

    if not os.path.exists(os.path.join(os.getcwd(), "luxetenebrae/runs", mode)):
        os.makedirs(os.path.join(os.getcwd(), "luxetenebrae/runs", mode))
    os.environ["COMPAS_LOGS_OUTPUT_DIR_PATH"] = os.path.join(os.getcwd(), "luxetenebrae/runs", mode) # This is the default output directory for COMPAS logs

    class pythonProgramOptions:
        """
        A class to store and access COMPAS program options in python
        """

        def __init__(self, config_file=DEFAULT_CONFIG_FILE, grid_filename=None,
                    random_seed_filename='randomSeed.txt', output_directory=None):
            # Do './COMPAS --help' to see all options
            # -- Define variables

            # environment variable COMPAS_EXECUTABLE_PATH is used for docker runs
            # if COMPAS_EXECUTABLE_PATH is not set (== None) we assume this is an
            # interactive run with python3
            # if COMPAS_EXECUTABLE_PATH is set (!= None) we assume this is a run
            # inside a docker container - we have different directories inside a
            # docker container (src, obj, bin), and the COMPAS executable resides
            # in the bin directory (rather than the src directory)

            # Load yaml file with options
            with open(config_file) as file:
                # The FullLoader parameter handles the conversion from YAML
                # scalar values to Python the dictionary format
                config = yaml.load(file, Loader=yaml.FullLoader)

            self.booleanChoices = config['booleanChoices'] if config['booleanChoices'] else {}
            self.numericalChoices = config['numericalChoices'] if config['numericalChoices'] else {}
            self.stringChoices = config['stringChoices'] if config['stringChoices'] else {}
            self.listChoices = config['listChoices'] if config['listChoices'] else {}

            compas_root_dir = os.environ.get('COMPAS_ROOT_DIR', REPO_ROOT)
            if compas_root_dir is None:
                warnings.warn(
                    'COMPAS_ROOT_DIR environment variable not set. Setting '
                    f'`export COMPAS_ROOT_DIR={REPO_ROOT}`'
                )
                os.environ['COMPAS_ROOT_DIR'] = REPO_ROOT
            compas_exe = os.path.join(compas_root_dir, 'src/COMPAS')
            compas_executable_override = os.environ.get('COMPAS_EXECUTABLE_PATH', compas_exe)
            print('compas_executable_override', compas_executable_override)
            self.compas_executable = compas_executable_override

            # If random_seed_filename is specified, overwrite the random seed from the yaml file
            if os.path.isfile(random_seed_filename):
                self.numericalChoices['--random-seed'] = int(np.loadtxt(random_seed_filename))

            # If grid is specified in pythonProgramOptions(), ignore the values from yaml file
            if grid_filename:
                self.grid_filename = grid_filename
                self.stringChoices['--grid'] = self.grid_filename
            else:
                if not self.stringChoices:
                    self.grid_filename = None
                else:
                    self.grid_filename = None
                    if '--grid' in self.stringChoices:
                        self.grid_filename = self.stringChoices['--grid']

            print('grid_filename', self.grid_filename)

            # environment variable COMPAS_LOGS_OUTPUT_DIR_PATH is used primarily for docker runs
            # if COMPAS_LOGS_OUTPUT_DIR_PATH is set (!= None) it is used as the value for the
            # --output-path option
            # if COMPAS_LOGS_OUTPUT_DIR_PATH is not set (== None) the current working directory
            # is used as the value for the --output-path option
            compas_logs_output_override = os.environ.get('COMPAS_LOGS_OUTPUT_DIR_PATH')

            if (compas_logs_output_override is None):
                self.stringChoices['--output-path'] = os.getcwd()
                # names the directory to be created and in which log files are created.  Default in COMPAS is "COMPAS_Output"
                self.stringChoices['--output-container'] = output_directory  
            else:
                self.stringChoices['--output-path'] = compas_logs_output_override
                self.stringChoices['--output-container'] = output_directory

                # environment variable COMPAS_INPUT_DIR_PATH is used primarily for docker runs
            # if COMPAS_INPUT_DIR_PATH is set (!= None) it is prepended to input filenames
            # (such as grid_filename and logfile_definitions)
            # if COMPAS_INPUT_DIR_PATH is not set (== None) the current working directory
            # is prepended to input filenames
            compas_input_path_override = os.environ.get('COMPAS_INPUT_DIR_PATH')

            if self.grid_filename != None:
                if compas_input_path_override == None:
                    self.grid_filename = os.getcwd() + '/' + self.grid_filename
                else:
                    self.grid_filename = compas_input_path_override + '/' + self.grid_filename

            if '--logfile-definitions' in self.stringChoices:
                self.logfile_definitions = self.stringChoices['--logfile-definitions']  # logfile record definitions file name (e.g. 'logdefs.txt')
            else:
                self.logfile_definitions = None

            if self.logfile_definitions != None:
                if compas_input_path_override == None:
                    self.logfile_definitions = os.getcwd() + '/' + self.logfile_definitions
                else:
                    self.logfile_definitions = compas_input_path_override + '/' + self.logfile_definitions

            self.makeCommandString()

        def makeCommandString(self):
            """
            This function generates a dictionary mapping COMPAS options to their specified
            values (or empty strings for boolean options). These are then combined into a string
            that can be run directly as a terminal command, or passed to the stroopwafel interface
            where some of them may be overwritten. Options not to be included in the command
            line should be set to pythons None (except booleans, which should be set to False)
            """

            ### Collect all options into a dictionary mapping option name to option value
            self.command = {'compas_executable': self.compas_executable}

            for boolKey, boolVal in self.booleanChoices.items():
                if boolVal is True:
                    self.command.update({boolKey: ''})
                elif boolVal is False:
                    self.command.update({boolKey: 'False'})

            for numKey, numVal in self.numericalChoices.items():
                if not numVal == None:
                    self.command.update({numKey: str(numVal)})

            for strKey, strVal in self.stringChoices.items():
                if not strVal == None:
                    self.command.update({strKey: strVal})

            for listKey, listVal in self.listChoices.items():
                if listVal:
                    self.command.update({listKey: ' '.join(map(str, listVal))})

            # Ensure the Compas executable is first, and not repeated.
            # Options are non-ordered.
            self.shellCommand = self.command['compas_executable']
            del self.command['compas_executable']
            for key, val in self.command.items():
                self.shellCommand += f' {key} {val}'

            return


    def runSubmit(cli_args=None, execute=True):
        parser = argparse.ArgumentParser(
            description='Run COMPAS using a config yaml (for settings refer to ./COMPAS --help)'
        )
        parser.add_argument('config_file', type=str, nargs="?", default=DEFAULT_CONFIG_FILE)
        parser.add_argument('--grid', type=str, default=None)
        args = parser.parse_args(cli_args)
        # -- Get the program options
        myoptions = pythonProgramOptions(config_file=args.config_file, grid_filename=args.grid)
        print(myoptions.shellCommand)
        if execute:  # Execute COMPAS shell string
            call(myoptions.shellCommand, shell=True)


    def main():
        cli_args = sys.argv[1:]
        runSubmit(cli_args=cli_args, execute=True)


    if __name__ == "__main__":
        main()
