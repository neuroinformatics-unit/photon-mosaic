import argparse
import shutil
from pathlib import Path

from calcium_imaging_automation.core.reader import ReadAquiredData
from calcium_imaging_automation.core.writer import DatashuttleWrapper


def setup(raw_data_path, folder_read_pattern, file_read_pattern, output_path):
    try:
        shutil.rmtree("/ceph/margrie/laura/cimaut/derivatives/")
        shutil.rmtree("/ceph/margrie/laura/cimaut/submitit/")
    except FileNotFoundError:
        print("No derivatives folder found")

    print(f"Reading data from {raw_data_path}")
    reader = ReadAquiredData(
        raw_data_path,
        folder_read_pattern,
        file_read_pattern,
    )
    print(f"Found {len(reader.datasets_paths)} datasets.")

    number_of_tiffs = reader.max_session_number(filetype="tif")
    print(f"Max of tiffs found: {number_of_tiffs}")

    writer = DatashuttleWrapper(output_path)
    print(f"Dataset names: {reader.dataset_names}")
    writer.create_folders(reader.dataset_names, session_number=number_of_tiffs)
    print("Folders created")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Example usage of the pipeline manager."
    )

    parser.add_argument(
        "raw_data_path", type=Path, help="Path to the raw data."
    )
    parser.add_argument(
        "output_path", type=Path, help="Path to the output data."
    )
    parser.add_argument(
        "--folder_read_pattern",
        type=str,
        help="Glob pattern for reading folder.",
        default="*",
    )
    parser.add_argument(
        "--file_read_pattern",
        type=str,
        help="List of glob patterns for reading files.",
        action="append",
    )

    args = parser.parse_args()

    try:
        setup(
            args.raw_data_path,
            args.folder_read_pattern,
            args.file_read_pattern,
            args.output_path,
        )

        print("Success")
    except Exception as e:
        print(f"Error: {e.args}")
        print(e.with_traceback(e.__traceback__))
