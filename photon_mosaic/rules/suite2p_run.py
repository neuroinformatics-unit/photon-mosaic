import datetime
import traceback
from pathlib import Path

import numpy as np
from snakemake.script import snakemake
from suite2p import run_s2p
from photon_mosaic.rules.s2p_options import get_edited_options

# Retrieve parameters and inputs from Snakemake
input_path = Path(snakemake.input[0])
ops_file = snakemake.input[1]
dataset_folder = Path(input_path).parent.parent

ops = get_edited_options(input_path, ops_file, dataset_folder)

db = {"data_path": input_path}
try:
    #  run suite2p
    ops_end = run_s2p(ops=ops)

    #  get registration metrics from ops
    metrics = {
        "regDX": ops_end.get("regDX", "NaN"),
        "regPC": ops_end.get("regPC", "NaN"),
        "tPC": ops_end.get("tPC", "NaN"),
    }

    #  append in the metrics file the new metrics
    with open(dataset_folder / "suite2p_metrics.txt", "w") as f:
        f.write("registration metrics: \n")
        for key, value in metrics.items():
            f.write(f"{key}: {value}\n")
    # make empty error file
    with open(dataset_folder / "error.txt", "a") as f:
        f.write("")
except Exception:
    with open(dataset_folder / "error.txt", "a") as f:
        #  add timestamp to the error file
        f.write(f"Error at {datetime.datetime.now()}\n")
        f.write(traceback.format_exc())
    with open(dataset_folder / "suite2p_metrics.txt", "w") as f:
        f.write("registration metrics: NaN\n")