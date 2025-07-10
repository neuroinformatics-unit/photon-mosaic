"""
Snakemake rule for running Poissson-numcodecs.
"""

from pathlib import Path
from typing import Optional

import numpy as np
from poisson_numcodecs import calibrate
import tifffile as tif
from pathlib import Path

from photon_mosaic.io.scanimage_utils import get_si_metadata, reshape_si_stack

def compress_movie(
        output_path: Path, 
        input_path: Path, 
        user_ops_dict: Optional[dict] = None,
    ):
    """
    This function converts a two-photon movie in a given dataset folder and saves the
    results in the specified paths. It also handles any exceptions
    that may occur during the process and logs them in an error
    file.

    Parameters
    ----------
    output_path : str
        The path to the folder where the compressed movie will be saved.
    input_path : Path
        The path to the TIFF file.
    user_ops_dict : dict, optional
        A dictionary containing user-provided options for Poission-numcodecs.

    Returns
    -------
    None
        The function compresses a two-photon movie and saves results to the specified paths.
    """
    if not user_ops_dict['numcodecs_ops']['scanimage']:
        raise NotImplementedError
    si_metadata = get_si_metadata(input_path)
    orig_file_name = input_path.name.split('.')[0]

    print("Reading tif image...")
    scan = tif.imread(input_path)

    if np.min(scan) < 0: # Make all values positive
        scan = scan + np.min(scan) * -1

    reshaped_scan = reshape_si_stack(scan, si_metadata)

    for plane in range(si_metadata['fpv']):
        plane_output_path = output_path / 'photon_flux' / 'plane{}'.format(plane)
        plane_output_path.mkdir(parents=True, exist_ok=True)

        for channel in range(si_metadata['nCh']):
            print("Converting movie to photon flux movie...")
            calibrator = calibrate.SequentialCalibratePhotons(reshaped_scan[:,plane,channel,:,:])

            if user_ops_dict['numcodecs_ops']['crop'] >0:
                calibrator.subsample_and_crop_video(crop=user_ops_dict['numcodecs_ops']['crop'])

            print("Calibrating photon sensitivity...")
            [photon_sensitivity, dark_signal] = calibrator.get_photon_sensitivity_parameters()
            print(f"Quantal size: {photon_sensitivity}\nIntercept: {dark_signal}\n")

            print("Getting photon flux movie...")
            print(
                "Performing the calculation photon flux movie = (scan  - dark_signal) / photon_sensitivity"
            )
            photon_counts_per_pixel_per_frame = calibrator.get_photon_flux_movie()

            saving_file_path = plane_output_path / (orig_file_name + '_chan{}.tif'.format(channel))
            tif.imwrite(saving_file_path, photon_counts_per_pixel_per_frame) 
            print(f"Saved photon flux movie to {str(saving_file_path)}")