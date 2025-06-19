"""
Contrast enhancement preprocessing step.

This module provides functions to enhance the contrast of images.
"""

from pathlib import Path
from typing import List, Optional, Union

import numpy as np
import tifffile
from skimage import exposure


def run(
    **kwargs,
) -> None:
    """
    Enhance the contrast of an image.

    Parameters
    ----------
    dataset_folder : Union[str, Path]
        Path to the dataset folder containing the input TIFF files.
    output_folder : Union[str, Path]
        Path to the output folder where the processed TIFF files will be saved.
    ses_idx : int
        Session index to process.
    glob_naming_pattern_tif : List[str]
        List of glob patterns to match TIFF files.
    output_path : Optional[Union[str, Path]], optional
        Path to save the enhanced image. If not provided, will be derived from
        output_folder.
    **kwargs : dict
        Additional keyword arguments:
        - percentile_low : float, optional
            Lower percentile for contrast stretching. Default is 1.
        - percentile_high : float, optional
            Upper percentile for contrast stretching. Default is 99.

    Returns
    -------
    None
        The function saves the enhanced image to the output path and returns
        nothing.
    """
    # Convert paths to Path objects
    dataset_folder = Path(kwargs["dataset_folder"])
    output_folder = Path(kwargs["output_folder"])
    tiff_file = dataset_folder / kwargs["tiff_name"]

    # Load the image
    img = tifffile.imread(tiff_file)

    # Get contrast parameters
    percentile_low = kwargs.get("percentile_low", 1)
    percentile_high = kwargs.get("percentile_high", 99)

    # Enhance contrast
    p_low, p_high = np.percentile(img, (percentile_low, percentile_high))
    img_enhanced = exposure.rescale_intensity(img, in_range=(p_low, p_high))

    # Append filename to output path
    output_path = output_folder / f"enhanced_{kwargs['tiff_name']}"

    # Save the enhanced image
    tifffile.imwrite(output_path, img_enhanced)
