import argparse
import datetime
import logging
from pathlib import Path

from calcium_imaging_automation.core.reader import ReadAllPathsInFolder
from calcium_imaging_automation.core.writer import DatashuttleWrapper


def main(
    raw_data_path: Path,
    output_path: Path,
    filetypes_of_interest: list,
    folder_read_pattern: str,
):
    """
    Draft usage of the pipeline, now consisting of read and write operations.
    """
    (output_path / "logs").mkdir(exist_ok=True)
    logging.basicConfig(
        # save also time anda date
        filename=str(
            output_path
            / "logs"
            / f"{datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')}.log"
        ),
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
    )

    reader = ReadAllPathsInFolder(
        raw_data_path, filetypes_of_interest, folder_read_pattern
    )
    logging.info(f"Found {len(reader.datasets_paths)} datasets.")
    logging.info(f"Dataset names: {reader.dataset_names}")

    writer = DatashuttleWrapper(output_path)

    number_of_tiffs = reader.max_session_number(filetype="tif")
    writer.create_folders(reader.dataset_names, session_number=number_of_tiffs)

    # [Placeholder for data processing]

    logging.info("Pipeline finished.")


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
    parser.add_argument(
        "--folder_read_pattern",
        type=str,
        help="Glob pattern for reading files.",
        default="*",
    )

    args = parser.parse_args()
    raw_data_path = args.raw_data_path
    output_path = args.output_path
    file_types = args.filetypes
    folder_read_pattern = args.folder_read_pattern

    main(raw_data_path, output_path, file_types, folder_read_pattern)
