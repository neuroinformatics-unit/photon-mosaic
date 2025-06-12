"""
Derotation preprocessing step for photon-mosaic.

This module provides a function to derotate image data using the derotation
package.
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
from derotation.derotate_batch import derotate

from .registry import register_step


@register_step("derotation")
def run(
    data: np.ndarray,
    dataset_folder: Optional[Path] = None,
    output_folder: Optional[Path] = None,
    **kwargs,
):
    """
    Derotate image data using the derotation pipeline.

    Parameters
    ----------
    data : np.ndarray
        The image data to derotate. Shape should be (time, height, width).
    dataset_folder : Path
        The path to the dataset folder.
    output_folder : Path
        The path to the output folder.
    **kwargs : dict
        Additional arguments for the derotation pipeline.

    Raises
    ------
    ValueError
        If required parameters are missing.
    """
    if dataset_folder is None or output_folder is None:
        raise ValueError(
            "dataset_folder and output_folder are required parameters"
        )

    # Run the derotation pipeline
    try:
        logging.info("Running derotation pipeline on image stack")
        derotate(
            dataset_folder=dataset_folder,
            output_folder=output_folder,
            **kwargs,
        )

    except Exception as e:
        logging.error(f"Derotation pipeline failed: {e}")
        raise
