from pathlib import Path
from typing import Dict, List

import numpy as np
from datashuttle.configs.config_class import Configs
from datashuttle.utils import folders
from PIL import Image


class DatashuttleWrapper:
    def __init__(self, output_path: Path) -> None:
        # This is supposed to run in the cluster and have direct access
        # to the central storages
        self.output_path = output_path
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
        # all_paths is a dictionary with keys: sub, ses
        self.all_paths: Dict[str, List[Path]] = folders.create_folder_trees(
            cfg=self.datashuttle_cfg,
            top_level_folder="derivatives",
            sub_names=[
                f"sub-{i}_{dataset_name}"
                for i, dataset_name in enumerate(dataset_names)
            ],
            ses_names=[f"ses-{i}" for i in range(session_number)],
            datatype="funcimg",
        )

    def get_dataset_path(self, dataset_name: str) -> Path:
        return next(
            (self.output_path / "derivatives").glob(f"*{dataset_name}*")
        )

    def save_image(
        self,
        image: np.ndarray,
        dataset_name: str,
        session_number: int,
        filename: str,
    ) -> Path:
        path = self.get_dataset_path(dataset_name)
        image = Image.fromarray(image).convert("L")
        image_path = (
            path / f"ses-{session_number}" / f"{filename}.png"
        )
        image.save(
            image_path,
            mode="PNG",
        )

        return image_path
