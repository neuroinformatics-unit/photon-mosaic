"""
Snakemake rule for preprocessing image data.

This module provides a function to run preprocessing on image data using the
preprocessing registry.
"""

from photon_mosaic.preprocessing import get_step


def run_preprocessing(
    output_path, config, dataset_folder=None, ses_idx=0, tiff_name=None
):
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

        kwargs["dataset_folder"] = dataset_folder
        kwargs["output_folder"] = output_path
        kwargs["ses_idx"] = ses_idx  # maybe not needed

        if step_name == "noop":
            kwargs["tiff_name"] = tiff_name

        # Add output path for contrast enhancement
        if step_name == "contrast":
            kwargs["output_path"] = output_path

        # Get the preprocessing function from the registry
        func = get_step(step_name)

        # Apply the preprocessing step
        func(**kwargs)
