#  Here you can explicitely change suite2p options
import numpy as np
from pathlib import Path

def get_edited_options(input_path, ops_file, dataset_folder):
        
    # load ops
    ops = np.load(ops_file, allow_pickle=True).item()

    #  Add the location in which the suite2p output will be saved
    ops["save_folder"] = str(dataset_folder)
    ops["save_path0"] = str(dataset_folder)
    ops["fast_disk"] = str(dataset_folder)
    ops["data_path"] = [input_path]

    #  change ops for non-rigid registration
    ops["nonrigid"] = True
    ops["block_size"] = [64, 64]
    ops["snr_thresh"] = 1.7
    ops["maxregshiftNR"] = 15