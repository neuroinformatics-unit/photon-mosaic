"""
Snakemake rule for preprocessing image data.

This module provides a function to run preprocessing on image data using the
preprocessing registry.
"""

from typing import Dict

from photon_mosaic.preprocessing import get_step


def get_input_files(dataset_folder, config, ses_idx):
    """
    Get input TIFF files based on the glob patterns defined in the config.

    Parameters
    ----------
    dataset_folder : Path
        Path to the dataset folder.
    config : dict
        Configuration dictionary containing glob patterns.

    Returns
    -------
    list of Path
        List of input TIFF file paths matching the patterns.
    """
    tiff_files = []

    # Iterate over all preprocessing steps
    for step in config["preprocessing"]["steps"]:
        # Get the glob pattern for this step
        glob_patterns = step["kwargs"].get("glob_naming_pattern_tif", [])

        # Add matching files for each pattern
        for pattern in glob_patterns:
            matched_files = [
                f for f in dataset_folder.rglob(pattern) if f.is_file()
            ]
            tiff_files.extend(matched_files)
            if not matched_files:
                raise FileNotFoundError(
                    f"No files found for pattern {pattern} in {dataset_folder}"
                )

    return tiff_files[ses_idx]


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


def run_preprocessing(output_path, config, dataset_folder=None, ses_idx=0):
    """
    Run preprocessing on image data.

    Parameters
    ----------
    output_path : str
        Path to save the preprocessed image.
    config : dict
        Configuration dictionary containing preprocessing steps.
    dataset_folder : str, optional
        Path to the dataset folder. This is needed for some preprocessing steps
        that require access to the dataset folder.
    ses_idx : int, optional
        Session index to process. Default is 0.

    Returns
    -------
    None
        The function saves the preprocessed data to the output path and returns
        nothing.
    """

    # Apply preprocessing steps
    for step in config["steps"]:
        step_name = step["name"]
        kwargs = step.get("kwargs", {})

        # Add dataset folder to kwargs if provided
        if dataset_folder and "dataset_folder" not in kwargs:
            kwargs["dataset_folder"] = dataset_folder
            kwargs["output_folder"] = output_path
            kwargs["ses_idx"] = ses_idx

        # Add output path for contrast enhancement
        if step_name == "contrast":
            kwargs["output_path"] = output_path

        # Get the preprocessing function from the registry
        func = get_step(step_name)

        # Apply the preprocessing step
        func(**kwargs)
