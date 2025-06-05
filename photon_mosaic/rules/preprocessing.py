"""
Snakemake rule for preprocessing image data.

This module provides a function to run preprocessing on image data using the
preprocessing registry.
"""

from typing import Dict

from tifffile import imread

from photon_mosaic.preprocessing import get_step


def get_output_pattern(tiff_name: str, config: Dict, tiff_paths: Dict) -> str:
    """
    Get the output pattern for a tiff file.

    Parameters
    ----------
    tiff_name : str
        The name of the tiff file (without extension)
    config : dict
        The configuration dictionary
    tiff_paths : dict
        Dictionary mapping dataset names to tiff file paths

    Returns
    -------
    str
        The output pattern with the tiff name substituted
    """
    # Check if we have a list of patterns
    if "output_patterns" in config["preprocessing"]:
        # Get the index of the current tiff in the list of tiffs
        tiff_index = list(tiff_paths.values()).index(tiff_name)
        # Get the corresponding pattern
        pattern = config["preprocessing"]["output_patterns"][tiff_index]
    else:
        # Use the single pattern
        pattern = config["preprocessing"]["output_pattern"]

    # Substitute the tiff name
    return pattern.format(tiff_name=tiff_name)


def run_preprocessing(input_path, output_path, config, dataset_folder=None):
    """
    Run preprocessing on image data.

    Parameters
    ----------
    input_path : Path
        Path to the input image file.
    output_path : str
        Path to save the preprocessed image.
    config : dict
        Configuration dictionary containing preprocessing steps.
    dataset_folder : str, optional
        Path to the dataset folder. This is needed for some preprocessing steps
        that require access to the dataset folder.
    """
    # Load images
    data = imread(input_path)

    # Apply preprocessing steps
    if "preprocessing" in config and "steps" in config["preprocessing"]:
        for step in config["preprocessing"]["steps"]:
            step_name = step["name"]
            kwargs = step.get("kwargs", {})

            # Add dataset folder to kwargs if provided
            if dataset_folder and "dataset_folder" not in kwargs:
                kwargs["dataset_folder"] = dataset_folder

            # Add output path for contrast enhancement
            if step_name == "contrast":
                kwargs["output_path"] = output_path

            # Get the preprocessing function from the registry
            func = get_step(step_name)

            # Apply the preprocessing step
            func(data, **kwargs)
