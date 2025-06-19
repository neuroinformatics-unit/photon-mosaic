"""
No-operation preprocessing step for photon-mosaic.

This module provides a function that returns the input data unchanged.
This is useful when preprocessing should be skipped.
"""

import shutil
from pathlib import Path


def run(
    dataset_folder: Path,
    output_folder: Path,
    tiff_name: str,
    **kwargs,
):
    """
    No-operation preprocessing step.

    Parameters
    ----------
    dataset_folder : Path
        Path to the dataset folder containing the input TIFF files.
    output_folder : Path
        Path to the output folder where the files will be copied.
    tiff_name : str
        Name of the TIFF file to copy.
    **kwargs : dict
        Additional keyword arguments (unused).

    Returns
    -------
    None
        The function copies the input TIFF file to the output directory
        without any modification and returns nothing.

    Notes
    -----
    The function will search for the TIFF file using rglob if it's not found
    at the expected location.
    """
    # Convert paths to Path objects if they're strings
    if isinstance(dataset_folder, str):
        dataset_folder = Path(dataset_folder)
    if isinstance(output_folder, str):
        output_folder = Path(output_folder)

    input_file = dataset_folder / tiff_name

    # Create output directory and copy file
    output_folder.mkdir(parents=True, exist_ok=True)
    try:
        shutil.copy2(input_file, output_folder / input_file.name)
    except FileNotFoundError:
        #  use rglob to find the correct path
        correct_path = next(dataset_folder.rglob(tiff_name))
        shutil.copy2(correct_path, output_folder / correct_path.name)
