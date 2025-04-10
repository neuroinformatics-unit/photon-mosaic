import traceback
from pathlib import Path

from suite2p import run_s2p

from photon_mosaic.s2p_options import get_edited_options


def run_suite2p(tiff_file, stat_path, bin_path, dataset_folder, user_ops_dict):
    save_folder = Path(stat_path).parent
    dataset_folder = Path(dataset_folder)
    ops = get_edited_options(
        input_path=dataset_folder,
        save_folder=save_folder,
        user_ops_dict=user_ops_dict,
    )
    try:
        ops_end = run_s2p(ops=ops)

        # save metrics, as before
        with open(dataset_folder / "suite2p_metrics.txt", "w") as f:
            f.write("registration metrics:\n")
            for key in ["regDX", "regPC", "tPC"]:
                f.write(f"{key}: {ops_end.get(key, 'NaN')}\n")

        # ensure output placeholders exist
        Path(stat_path).touch()
        Path(bin_path).touch()

    except Exception as e:
        with open(dataset_folder / "error.txt", "a") as f:
            f.write(f"Error: {e}\n")
            f.write(traceback.format_exc())
