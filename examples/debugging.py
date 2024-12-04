import shutil
from pathlib import Path

from derotation.analysis.metrics import stability_of_most_detected_blob
from derotation.derotate_batch import derotate

from calcium_imaging_automation.core.pipeline import orchestrator

try:
    shutil.rmtree("/ceph/margrie/laura/cimaut/derivatives/")
    shutil.rmtree("/ceph/margrie/laura/cimaut/submitit/")
except FileNotFoundError:
    print("No derivatives folder found")

orchestrator(
    raw_data_path=Path(
        "/nfs/winstor/margrie/SimonWeiler/RawData/Invivo_imaging/3photon_rotation/shared/"
    ),
    output_path=Path("/ceph/margrie/laura/cimaut/"),
    folder_read_pattern="2*",
    file_read_pattern=["rotation_00001.tif", "*.bin"],
    experiment_name="submitit_04",
    preprocessing_function=derotate,
    compute_metric=stability_of_most_detected_blob,
    # suite2p_ops_path="/ceph/margrie/laura/derotation/suite2p/laura_ops.npy",
)
