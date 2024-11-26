from pathlib import Path
from typing import List


class ReadAquiredData:
    def __init__(
        self,
        raw_data_folder: Path,
        folder_read_pattern: str,
        file_read_pattern: List[str],
    ):
        """
        Class to handle filepaths and dataset names in the raw data folder.
        It can load folders and files based on the provided patterns, allowing
        flexibility in the data structure of origin.
        It also provides the maximum number of sessions for each dataset based
        on the total number of files found in the dataset folders, by default
        it searches for tif files.

        Parameters
        ----------
        raw_data_folder : Path
            The path to the raw data folder.
        folder_read_pattern : str
            The pattern to search for folders in the raw data folder. It
            corresponds to the naming convention of the datasets.
        file_read_pattern : List[str]
            The patterns to search for files in the dataset folders. It
            corresponds to the naming convention of the files in the dataset
            folders.
        """
        self.folder_read_pattern = folder_read_pattern
        self.file_read_pattern = file_read_pattern

        self.datasets_paths = self.get_folders_first_layer(raw_data_folder)
        self.dataset_names = [
            dataset_path.name for dataset_path in self.datasets_paths
        ]

    def get_folders_first_layer(self, file_path: Path) -> List[Path]:
        """
        Get the first layer of folders in the raw data folder. The rest
        of the class assumes that the first layer of folders corresponds
        to the dataset folders.

        Parameters
        ----------
        file_path : Path
            The path to the raw data folder.

        Returns
        -------
        List[Path]
            The list of paths to the dataset folders.
        """
        return list(file_path.glob(self.folder_read_pattern))

    def get_files_paths_by_format(
        self, folder: Path, filetype="tif"
    ) -> List[Path]:
        """
        Get the paths to the files in the dataset folders based on the
        provided file type. By default, it searches for tif files.

        Parameters
        ----------
        folder : Path
            The path to the dataset folder.
        filetype : str, optional
            The file type to search for in the dataset folder, by default
            "tif".

        Returns
        -------
        List[Path]
            The list of paths to the files in the dataset folder.
        """
        return list(folder.rglob(filetype))

    def total_objects_by_extension(self, folder: Path) -> dict:
        """
        Get the total number of files in the dataset folder based on the
        extensions included in the file_read_pattern.

        Parameters
        ----------
        folder : Path
            The path to the dataset folder.

        Returns
        -------
        dict
            The dictionary with the number of files for each extension in the
            patterns found in file_read_pattern.
        """

        return {
            filetype.split(".")[-1]: len(
                self.get_files_paths_by_format(folder, filetype)
            )
            for filetype in self.file_read_pattern
        }

    def max_session_number(self, filetype="tif", max_allowed=1) -> int:
        """
        Get the maximum number of sessions for each dataset based on the total
        number of files found in the dataset folders. By default, it searches
        for tif files and allows a maximum of 5 sessions. It assumes that every
        tif file corresponds to an experimental session.

        Parameters
        ----------
        filetype : str, optional
            The file type to search for in the dataset folder, by default
            "tif".
        max_allowed : int, optional
            The maximum number of sessions allowed, by default 5.

        Returns
        -------
        int
            The maximum number of sessions for each dataset.
        """

        total_tif_number = [
            self.total_objects_by_extension(dataset_path).get(filetype, 0)
            for dataset_path in self.datasets_paths
        ]

        return min(max(total_tif_number), max_allowed)
