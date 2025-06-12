"""
No-operation preprocessing step for photon-mosaic.

This module provides a function that returns the input data unchanged.
This is useful when preprocessing should be skipped.
"""

import shutil
from pathlib import Path

from .registry import register_step


@register_step("noop")
def run(data=None, **kwargs):
    """
    No-operation preprocessing step.

    Parameters
    ----------
    data : np.ndarray or None, optional
        The image data. If None, the function will just create the output
        directory.
    **kwargs : dict
        Additional arguments. If data is None, kwargs should contain:
        - dataset_folder: Path to the input file
        - output_folder: Path to save the output
        - ses_idx: Session index

    Returns
    -------
    None
        The function either returns the input data unchanged or copies the input file
        to the output directory and returns nothing.
    """
    if data is None:
        # If no data provided, just copy the input file to output
        dataset_folder = Path(kwargs["dataset_folder"])
        output_folder = Path(kwargs["output_folder"])
        ses_idx = kwargs["ses_idx"]

        # Get the input file from the glob pattern
        glob_pattern = kwargs["glob_naming_pattern_tif"][ses_idx]
        input_file = next(dataset_folder.rglob(glob_pattern))

        # Create output directory and copy file
        output_folder.mkdir(parents=True, exist_ok=True)
        shutil.copy2(input_file, output_folder / input_file.name)
