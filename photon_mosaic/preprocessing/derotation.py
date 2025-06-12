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
    ses_idx: int = 0,
    **kwargs,
):
    """
    Derotate image data using the derotation pipeline.

    Parameters
    ----------
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
    ses_idx : int, optional
        Session index to process. Default is 0.
    **kwargs : dict
        Additional arguments for the derotation pipeline.

    Returns
    -------
    None
        The function saves the derotated data to the output folder and returns nothing.

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
        pattern_tif = glob_naming_pattern_tif
    else:
        pattern_tif = glob_naming_pattern_tif[ses_idx]

    if isinstance(glob_naming_pattern_bin, str):
        pattern_bin = glob_naming_pattern_bin
    else:
        pattern_bin = glob_naming_pattern_bin[ses_idx]

    logging.info(
        f"Running derotation pipeline for {pattern_tif} with" f" {pattern_bin}"
    )
    derotate(
        dataset_folder=dataset_folder,
        output_folder=output_folder,
        glob_naming_pattern_tif=pattern_tif,
        glob_naming_pattern_bin=pattern_bin,
        folder_suffix="incremental" if "increment" in pattern_tif else "full",
        **kwargs,
    )
