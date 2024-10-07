from pathlib import Path
from typing import List

from datashuttle.configs.config_class import Configs
from datashuttle.utils import folders


class DatashuttleWrapper:
    def __init__(self, output_path: Path) -> None:
        # This is supposed to run in the cluster and have direct access
        # to the central storages
        self.datashuttle_cfg = Configs(
            project_name=output_path.name,
            file_path=output_path,
            input_dict={
                "local_path": output_path,
                "central_path": "",
                "connection_method": "local_filesystem",
            },
        )

    def create_folders(self, dataset_names: List[str], session_number) -> None:
        folders.create_folder_trees(
            cfg=self.datashuttle_cfg,
            top_level_folder="derivatives",
            sub_names=[
                f"sub-{i}_{dataset_name}"
                for i, dataset_name in enumerate(dataset_names)
            ],
            ses_names=[f"ses-{i}" for i in range(session_number)],
            datatype="funcimg",
        )
