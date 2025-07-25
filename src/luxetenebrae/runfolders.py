import os
from datetime import datetime


def get_run_folder(base_dir="Runs", run_id=None):
    """
    Create or get a run folder. If run_id is None, use current timestamp.
    Updates latest.txt pointer in base_dir.
    Returns the absolute path to the run folder.
    """
    os.makedirs(base_dir, exist_ok=True)
    if run_id is None:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_folder = os.path.join(base_dir, f"run_{run_id}")
    os.makedirs(run_folder, exist_ok=True)
    # Update latest.txt pointer
    with open(os.path.join(base_dir, "latest.txt"), "w") as f:
        f.write(os.path.abspath(run_folder))
    return os.path.abspath(run_folder)


def get_latest_run_folder(base_dir="Runs"):
    """
    Return the absolute path to the latest run folder (from latest.txt).
    Raises RuntimeError if not found.
    """
    latest_path = os.path.join(base_dir, "latest.txt")
    if not os.path.exists(latest_path):
        raise RuntimeError("No previous run found. Please run the first stage.")
    with open(latest_path) as f:
        return f.read().strip()
