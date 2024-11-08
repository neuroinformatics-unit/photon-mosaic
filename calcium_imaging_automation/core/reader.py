from pathlib import Path
from typing import List


class ReadAllPathsInFolder:
    def __init__(
        self,
        raw_data_folder: Path,
        folder_read_pattern: str,
        file_read_pattern: list,
    ):
        self.folder_read_pattern = folder_read_pattern
        self.file_read_pattern = file_read_pattern

        self.datasets_paths = self.get_folders_first_layer(raw_data_folder)
        self.dataset_names = [
            dataset_path.name for dataset_path in self.datasets_paths
        ]

    def get_folders_first_layer(self, file_path: Path) -> List[Path]:
        return list(file_path.glob(self.folder_read_pattern))

    def get_files_paths(self, folder: Path) -> List[Path]:
        return [
            file
            for pattern in self.folder_read_pattern
            for file in folder.rglob(pattern)
        ]

    def total_objects_by_format(self, folder: Path) -> dict:
        return {
            filetype: len(self.get_files_paths(folder))
            for filetype in self.folder_read_pattern.split(".")[-1]
        }

    def max_session_number(self, filetype="tif", max_allowed=1) -> int:
        total_tif_number = [
            self.total_objects_by_format(dataset_path).get(filetype, 0)
            for dataset_path in self.datasets_paths
        ]

        return min(max(total_tif_number), max_allowed)
