import argparse
from pathlib import Path

from calcium_imaging_automation.core.reader import ReadAllPathsInFolder
from calcium_imaging_automation.core.writer import DatashuttleWrapper


def main(raw_data_path: Path, output_path: Path, filetypes_of_interest: list):
    """
    Draft usage of the pipeline, now consisting of read and write operations.
    """
    reader = ReadAllPathsInFolder(raw_data_path, filetypes_of_interest)

    writer = DatashuttleWrapper(output_path)
    number_of_tiffs = reader.max_session_number(filetype="tif")
    writer.create_folders(reader.dataset_names, session_number=number_of_tiffs)

    # [Placeholder for data processing]


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
        "--filetypes",
        type=list,
        nargs="+",
        help="Filetypes of interest.",
        default=["tif", "bin"],
    )

    args = parser.parse_args()
    raw_data_path = args.raw_data_path
    output_path = args.output_path
    file_types = args.filetypes

    main(raw_data_path, output_path, file_types)
