"""
Derotation preprocessing step for photon-mosaic.

This module provides a function to derotate image data using the derotation
package.
"""

import logging
from pathlib import Path
from typing import List, Optional, Union

from derotation.derotate_batch import derotate

from .registry import register_step


@register_step("derotation")
def run(
    dataset_folder: Optional[Path] = None,
    output_folder: Optional[Path] = None,
    glob_naming_pattern_tif: Union[str, List[str]] = "*.tif",
    glob_naming_pattern_bin: Union[str, List[str]] = "*.bin",
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
    glob_naming_pattern_tif : Union[str, List[str]]
        Pattern(s) to match tif files. Can be a single pattern or a list of
        specific filenames.
    glob_naming_pattern_bin : Union[str, List[str]]
        Pattern(s) to match bin files. Can be a single pattern or a list of
        specific filenames.
    **kwargs : dict
        Additional arguments for the derotation pipeline.

    Raises
    ------
    ValueError
        If required parameters are missing or if tif/bin patterns don't match
        in length.
    """
    if dataset_folder is None or output_folder is None:
        raise ValueError(
            "dataset_folder and output_folder are required parameters"
        )

    # Convert string patterns to lists if needed
    if isinstance(glob_naming_pattern_tif, str):
        glob_naming_pattern_tif = [glob_naming_pattern_tif]
    if isinstance(glob_naming_pattern_bin, str):
        glob_naming_pattern_bin = [glob_naming_pattern_bin]

    # Ensure we have matching numbers of tif and bin patterns
    if len(glob_naming_pattern_tif) != len(glob_naming_pattern_bin):
        raise ValueError(
            f"Number of tif patterns ({len(glob_naming_pattern_tif)}) "
            "must match number of bin patterns "
            f"({len(glob_naming_pattern_bin)})"
        )

    # Run the derotation pipeline for each pair
    try:
        for tif_pattern, bin_pattern in zip(
            glob_naming_pattern_tif, glob_naming_pattern_bin
        ):
            logging.info(
                f"Running derotation pipeline for {tif_pattern} with"
                f" {bin_pattern}"
            )
            derotate(
                dataset_folder=dataset_folder,
                output_folder=output_folder,
                glob_naming_pattern_tif=tif_pattern,
                glob_naming_pattern_bin=bin_pattern,
                folder_suffix="incremental"
                if "increment" in tif_pattern
                else "full",
                **kwargs,
            )

    except Exception as e:
        logging.error(f"Derotation pipeline failed: {e}")
        raise
