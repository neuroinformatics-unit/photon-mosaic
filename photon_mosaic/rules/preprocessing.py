"""
Snakemake rule for preprocessing image data.

This module provides a function to run preprocessing on image data using the
preprocessing registry.
"""

import numpy as np
from tifffile import imread

from photon_mosaic.preprocessing import get_step


def run_preprocessing(input_paths, output_path, config, dataset_folder=None):
    """
    Run preprocessing on image data.

    Parameters
    ----------
    input_paths : list
        List of paths to input image files.
    output_path : str
        Path to save the preprocessed image.
    config : dict
        Configuration dictionary containing preprocessing steps.
    dataset_folder : str, optional
        Path to the dataset folder. This is needed for some preprocessing steps
        that require access to the dataset folder.
    """
    # Load images
    images = [imread(p) for p in input_paths]
    data = np.stack(images)

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
