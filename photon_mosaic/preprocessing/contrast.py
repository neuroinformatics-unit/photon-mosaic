"""
Contrast enhancement preprocessing step.

This module provides functions to enhance the contrast of images.
"""

from pathlib import Path
from typing import List, Optional, Union

import numpy as np
import tifffile
from skimage import exposure

from .registry import register_step


@register_step("contrast")
def run(
    dataset_folder: Union[str, Path],
    output_folder: Union[str, Path],
    ses_idx: int,
    glob_naming_pattern_tif: List[str],
    output_path: Optional[Union[str, Path]] = None,
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
    dataset_folder = Path(dataset_folder)
    output_folder = Path(output_folder)
    if output_path is not None:
        output_path = Path(output_path)

    # Find the input TIFF file
    tiff_files = []
    for pattern in glob_naming_pattern_tif:
        matched_files = [
            f for f in dataset_folder.rglob(pattern) if f.is_file()
        ]
        tiff_files.extend(matched_files)
        if not matched_files:
            raise FileNotFoundError(
                f"No files found for pattern {pattern} in {dataset_folder}"
            )

    # Get the specific session file
    input_path = tiff_files[ses_idx]

    # Load the image
    img = tifffile.imread(input_path)

    # Get contrast parameters
    percentile_low = kwargs.get("percentile_low", 1)
    percentile_high = kwargs.get("percentile_high", 99)

    # Enhance contrast
    p_low, p_high = np.percentile(img, (percentile_low, percentile_high))
    img_enhanced = exposure.rescale_intensity(img, in_range=(p_low, p_high))

    # Append filename to output path
    output_path = output_folder / f"enhanced_{input_path.name}"

    # Save the enhanced image
    tifffile.imwrite(output_path, img_enhanced)
