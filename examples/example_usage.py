import argparse
import datetime
import logging
from pathlib import Path
from typing import List

import mlflow
import numpy as np

from calcium_imaging_automation.core.reader import ReadAllPathsInFolder
from calcium_imaging_automation.core.writer import DatashuttleWrapper


def main(
    raw_data_path: Path,
    output_path: Path,
    folder_read_pattern: str,
    file_read_pattern: List[str],
):
    """
    Draft usage of the pipeline, now consisting of read and write operations.
    """
    # --- Setup experiment-wide logging to file ---
    (output_path / "logs").mkdir(exist_ok=True)
    logging.basicConfig(
        # Save also time and date
        filename=str(
            output_path
            / "logs"
            / f"{datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')}.log"
        ),
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
    )

    # --- Read folders and files ---
    reader = ReadAllPathsInFolder(
        raw_data_path,
        folder_read_pattern,
        file_read_pattern,
    )
    logging.info(f"Found {len(reader.datasets_paths)} datasets.")
    logging.info(f"Dataset names: {reader.dataset_names}")

    number_of_tiffs = reader.max_session_number(filetype="tif")
    logging.info(f"Max of tiffs found: {number_of_tiffs}")

    # --- Write folders and files ---
    writer = DatashuttleWrapper(output_path)
    writer.create_folders(reader.dataset_names, session_number=number_of_tiffs)

    for dataset in reader.datasets_paths:
        dataset_name = dataset.stem
        for session in range(1, number_of_tiffs + 1):
            with (
                mlflow.start_run()
            ):  # Start a new MLflow run for each dataset-session
                # Log session-specific parameters
                mlflow.log_param("dataset_name", dataset_name)
                mlflow.log_param("session_number", session)
                mlflow.log_param("raw_data_path", str(raw_data_path))
                mlflow.log_param("output_path", str(output_path))
                mlflow.log_param("folder_read_pattern", folder_read_pattern)
                mlflow.log_param("file_read_pattern", file_read_pattern)

                logging.info(
                    f"Processing dataset {dataset_name} session {session}..."
                )

                # Mock processing
                data = np.random.rand(100, 100)
                metric_measured = np.random.rand()

                # Log metric with MLflow
                mlflow.log_metric("metric_measured", metric_measured)

                # Save image in session folder
                writer.save_image(
                    image=data,
                    run_id=session,
                    dataset_name=dataset_name,
                    session_number=session,
                    filename="image",
                )

                # Log that the run is complete for this session
                logging.info(
                    f"Completed MLflow run for dataset {dataset_name} "
                    + f"session {session}"
                )

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
    raw_data_path = args.raw_data_path
    output_path = args.output_path
    folder_read_pattern = args.folder_read_pattern
    file_read_pattern = args.file_read_pattern

    main(
        raw_data_path,
        output_path,
        folder_read_pattern,
        file_read_pattern,
    )
