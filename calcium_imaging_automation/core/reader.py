from pathlib import Path
from typing import List


class ReadAllPathsInFolder:
    def __init__(
        self,
        raw_data_folder: Path,
        filetypes_of_interest: List[str],
        folder_read_pattern: str,
    ):
        self.filetypes_of_interest = filetypes_of_interest
        self.folder_read_pattern = folder_read_pattern

        self.datasets_paths = self.get_folders_first_layer(raw_data_folder)
        self.dataset_names = [
            dataset_path.name for dataset_path in self.datasets_paths
        ]

    def get_folders_first_layer(self, file_path: Path) -> List[Path]:
        return list(file_path.glob(self.folder_read_pattern))

    def get_files_paths(self, folder: Path) -> List[Path]:
        return [
            file
            for filetype in self.filetypes_of_interest
            for file in folder.rglob(f"*.{filetype}")
        ]

    def total_objects_by_filetype(self, folder: Path) -> dict:
        return {
            filetype: len(self.get_files_paths(folder))
            for filetype in self.filetypes_of_interest
        }

    def max_session_number(self, filetype="tif", max_allowed=5) -> int:
        total_tif_number = [
            self.total_objects_by_filetype(dataset_path).get(filetype, 0)
            for dataset_path in self.datasets_paths
        ]

        return min(max(total_tif_number), max_allowed)
