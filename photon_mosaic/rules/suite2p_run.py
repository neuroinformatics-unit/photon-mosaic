import traceback
from pathlib import Path
from typing import Optional

from suite2p import run_s2p

from photon_mosaic.s2p_options import get_edited_options


def run_suite2p(
    tiff_file: str,
    stat_path: str,
    bin_path: str,
    dataset_folder: Path,
    user_ops_dict: Optional[dict] = None,
):
    """
    This function runs Suite2P on a given TIFF file and saves the
    results in the specified paths. It also handles any exceptions
    that may occur during the process and logs them in an error
    file.

    Parameters
    ----------
    tiff_file : str
        The path to the TIFF file to be processed.
    stat_path : str
        The path where the Suite2P statistics will be saved.
    bin_path : str
        The path where the Suite2P binary data will be saved.
    dataset_folder : str
        The path to the folder containing the dataset.
    user_ops_dict : dict, optional
        A dictionary containing user-provided options to override
        the default Suite2P options. The default is None.
    """
    save_folder = Path(stat_path).parents[1]

    ops = get_edited_options(
        input_path=dataset_folder,
        save_folder=save_folder,
        user_ops_dict=user_ops_dict,
    )
    try:
        run_s2p(ops=ops)
    except Exception as e:
        with open(dataset_folder / "error.txt", "a") as f:
            f.write(f"Error: {e}\n")
            f.write(traceback.format_exc())
