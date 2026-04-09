import os, sys
import math
import yaml
import numpy as np
import warnings
import argparse
from subprocess import call
import secrets
from compas_python_utils.preprocessing.runSubmit import runSubmit
import os
os.environ["COMPAS_EXECUTABLE_PATH"] = "/data/a.saricaoglu/repo/luxetenebrae/COMPAS/src/bin/COMPAS"
os.environ["COMPAS_ROOT_DIR"] = "/data/a.saricaoglu/repo/luxetenebrae/COMPAS"
class GenerateBinaryPopulation:
    """
    Runs COMPAS to generate a population of binary systems based on input data.
    """

    def __init__(self, mode, number_of_systems, starting_seed, config_files=None, output_dir=None):
        self.mode = mode
        # One COMPAS run ~ 1000 systems
        self.number_of_runs = int(math.ceil(number_of_systems / 1000.0))
        self.starting_seed = starting_seed
        # Allow a single path or a list of paths
        if isinstance(config_files, str):
            config_files = [config_files]
        self.config_files = config_files
        self.output_dir = output_dir  # optional override

    def _prepare_config(self, config_path, seed, run_idx, out_dir):
        """
        Copy config file and override --random-seed with the given seed.
        Returns the path to the modified config file.
        """
        with open(config_path, 'r') as f:
            cfg = yaml.safe_load(f)

        # Set the random seed in numericalChoices
        if 'numericalChoices' not in cfg:
            cfg['numericalChoices'] = {}
        cfg['numericalChoices']['--random-seed'] = seed

        # Save modified config to output dir
        modified_cfg_path = os.path.join(out_dir, f"config_run{run_idx}_seed{seed}.yaml")
        with open(modified_cfg_path, 'w') as f:
            yaml.dump(cfg, f)

        print(f"[COMPAS] Prepared config: {modified_cfg_path} with seed={seed}")
        return modified_cfg_path
    
    def _default_config_for_mode(self):
        here = os.path.dirname(__file__)
        repo_root = os.path.abspath(os.path.join(here, "../../"))
        return os.path.join(repo_root, f"config/config{self.mode}.yaml"), repo_root
    
    def _write_random_seed(self, out_dir):
        """Create a new randomSeed.txt in out_dir and append to seeds_history.txt."""
        os.makedirs(out_dir, exist_ok=True)

        # If an existing seed file is present, archive its content
        seed_path = os.path.join(out_dir, "randomSeed.txt")
        history_path = os.path.join(out_dir, "seeds_history.txt")

        old_seed = None
        if os.path.exists(seed_path):
            try:
                old_seed = open(seed_path, "r").read().strip()
            except Exception:
                old_seed = None

        # Generate a new seed
        seed = secrets.randbits(32)

        # Append old seed and new seed to history (timestamped)
        import datetime as dt
        with open(history_path, "a") as hist:
            ts = dt.datetime.utcnow().isoformat(timespec="seconds") + "Z"
            if old_seed:
                hist.write(f"{ts} previous {old_seed}\n")
            hist.write(f"{ts} new {seed}\n")

        # Write the new seed to randomSeed.txt
        with open(seed_path, "w") as f:
            f.write(str(seed))

        print(f"[COMPAS] Wrote random seed {seed} to {seed_path} (history: {history_path})")
        return seed_path


    def _call_runSubmit(self, config_path):
        """
        Call runSubmit with CLI-like arguments. Falls back to sys.argv if needed.
        """
        try:
            # Most implementations accept cli_args=[...]
            return runSubmit(cli_args=[config_path])
        except TypeError:
            # Fallback: simulate CLI
            argv_bak = sys.argv[:]
            try:
                sys.argv = ["runSubmit.py", config_path]
                return runSubmit()
            finally:
                sys.argv = argv_bak
        print(f'arguments restored: {sys.argv}')

    def run_COMPAS_simulation(self):
        default_cfg, repo_root = self._default_config_for_mode()

        # Prepare output directory env
        out_dir = self.output_dir or os.path.join(repo_root, "runs", self.mode)
        os.makedirs(out_dir, exist_ok=True)
        os.environ["COMPAS_LOGS_OUTPUT_DIR_PATH"] = out_dir

        # Determine which configs to use
        cfgs = self.config_files or [default_cfg]
        seed = self.starting_seed
        for run_idx in range(self.number_of_runs):
            cfg = cfgs[run_idx % len(cfgs)]
            print(f"[COMPAS] Run {run_idx+1}/{self.number_of_runs} using config: {cfg}, seed: {seed}")
            if not os.path.isfile(cfg):
                raise FileNotFoundError(f"Config file not found: {cfg}")

            # Prepare a modified config with the seed baked in
            modified_cfg = self._prepare_config(cfg, seed, run_idx, out_dir)
            self._call_runSubmit(modified_cfg)
            seed += 1000

        # Save for introspection
        self._last_repo_root = repo_root
        self._last_default_config = default_cfg
        self._last_output_dir = out_dir

    def get_info(self):
        return {
            "mode": self.mode,
            "number_of_systems": self.number_of_runs * 1000,
            "repository_root": getattr(self, "_last_repo_root", None),
            "default_config_file": getattr(self, "_last_default_config", None),
            "output_directory": getattr(self, "_last_output_dir", os.environ.get("COMPAS_LOGS_OUTPUT_DIR_PATH")),
        }