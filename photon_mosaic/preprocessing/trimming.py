"""
Trimming preprocessing step.

This preprocessing step trims unwanted image pixels at the edge of the original files.
"""

import logging
from pathlib import Path

import tifffile

logger = logging.getLogger(__name__)


def run(
    dataset_folder: Path,
    output_folder: Path,
    tiff_name: str,
    trim_x: int = 5,
    trim_y: int = 5,
    **kwargs,
) -> None:
    """
    Trim the edges of an image.

    Parameters
    ----------
    dataset_folder : Path
        Path to the dataset folder containing the input TIFF files.
    output_folder : Path
        Path to the output folder where the processed TIFF files will be saved.
    tiff_name : str
        Name of the TIFF file to process.
    trim_x : int, optional
        Number of pixels to be trimmed for X axis. Default is 5.
    trim_y : int, optional
        Number of pixels to be trimmed for Y axis. Default is 5.
    **kwargs : dict
        Additional keyword arguments (unused).

    Returns
    -------
    None
        The function saves the trimmed image to the output folder with the
        prefix "trimmed_" and returns nothing.

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

    tiff_file = dataset_folder / tiff_name

    # Load the image
    try:
        img = tifffile.imread(tiff_file)
    except FileNotFoundError:
        #  use rglob to find the correct path
        correct_path = next(dataset_folder.rglob(tiff_name))
        img = tifffile.imread(correct_path)

    # Trim pixels
    if img.ndim == 3:
        img_trimmed = img[:, trim_x:-trim_x, trim_y:-trim_y]
    elif img.ndim == 4:
        img_trimmed = img[:, :, trim_x:-trim_x, trim_y:-trim_y]
    else:
        raise NotImplementedError(
            "Trimming preprocessing only supports 3D or 4D TIFF stacks"
        )

    # Append filename to output path
    output_path = output_folder / f"trimmed_{tiff_name}"

    # Save the enhanced image
    tifffile.imwrite(output_path, img_trimmed)
