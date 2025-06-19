"""
No-operation preprocessing step for photon-mosaic.

This module provides a function that returns the input data unchanged.
This is useful when preprocessing should be skipped.
"""

import shutil
from pathlib import Path


def run(**kwargs):
    """
    No-operation preprocessing step.

    Parameters
    ----------
    **kwargs : dict
        Additional arguments. If data is None, kwargs should contain:
        - dataset_folder: Path to the input file
        - output_folder: Path to save the output
        - ses_idx: Session index

    Returns
    -------
    None
        The function either returns the input data unchanged or copies the
        input file to the output directory and returns nothing.
    """
    # If no data provided, just copy the input file to output
    dataset_folder = Path(kwargs["dataset_folder"])
    output_folder = Path(kwargs["output_folder"])

    input_file = dataset_folder / kwargs["tiff_name"]

    # Create output directory and copy file
    output_folder.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_file, output_folder / input_file.name)
