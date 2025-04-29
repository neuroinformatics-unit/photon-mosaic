from pathlib import Path
from typing import Optional

from suite2p import default_ops


def get_edited_options(
    input_path: Path, save_folder: Path, user_ops_dict: Optional[dict] = None
) -> dict:
    """This function generates a dictionary of options for Suite2P
    by loading the default options and then modifying them with
    user-provided options. The function also sets the required
    runtime paths for saving the results.

    Parameters
    ----------
    input_path : Path
        The path to the input data folder.
    save_folder : Path
        The path to the folder where the results will be saved.
    user_ops_dict : dict, optional
        A dictionary containing user-provided options to override
        the default options. The default is None.

    Returns
    -------
    dict
        A dictionary containing the Suite2P options, including
        the user-provided options and the required runtime paths.

    Raises
    ------
    ValueError
        If a user-provided option is not valid for Suite2P.
    """

    ops = default_ops()

    # Override with user-provided subset of keys
    if user_ops_dict:
        for key, val in user_ops_dict.items():
            if key not in ops:
                raise ValueError(f"Invalid Suite2p option: {key}")
            ops[key] = val

    # Add required runtime paths
    ops["save_folder"] = str(save_folder)
    ops["save_path0"] = str(save_folder)
    ops["fast_disk"] = str(save_folder.parent)
    ops["data_path"] = [str(input_path)]

    return ops
