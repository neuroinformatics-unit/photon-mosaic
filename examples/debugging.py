from pathlib import Path

from calcium_imaging_automation.core.pipeline import pipeline

pipeline(
    raw_data_path=Path(
        "/nfs/winstor/margrie/SimonWeiler/RawData/Invivo_imaging/3photon_rotation/shared/"
    ),
    output_path=Path("/ceph/margrie/laura/cimaut/"),
    folder_read_pattern="2*",
    file_read_pattern=["rotation_00001.tif", "*.bin"],
)
