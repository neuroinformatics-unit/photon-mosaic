import argparse
import datetime
import logging
from pathlib import Path
from typing import List
import setuptools_scm

import mlflow
import numpy as np

from calcium_imaging_automation.core.reader import ReadAquiredData
from calcium_imaging_automation.core.writer import DatashuttleWrapper


def main(
    raw_data_path: Path,
    output_path: Path,
    folder_read_pattern: str,
    file_read_pattern: List[str],
    experiment_name: str = "pipeline_test",
):
    # --- Setup experiment-wide logging to file ---
    (output_path / "logs").mkdir(exist_ok=True)
    logging.basicConfig(
        filename=str(
            output_path
            / "logs"
            / f"{datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')}.log"
        ),
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
    )

    # --- Setup MLflow tracking ---
    mlflow_tracking_dir = output_path / "mlflow"
    mlflow.set_tracking_uri(str(mlflow_tracking_dir))
    mlflow.set_experiment(experiment_name)

    # --- Read folders and files ---
    reader = ReadAquiredData(
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
        for session in range(0, number_of_tiffs):
            # Generate mock data
            data = np.random.rand(100, 100)

            # Start a new MLflow experiment for each dataset-session
            with mlflow.start_run() as parent_run:
                # Log session-specific parameters
                mlflow.log_param("dataset_name", dataset_name)
                mlflow.log_param("session_number", session)
                mlflow.log_param("raw_data_path", str(raw_data_path))
                mlflow.log_param("output_path", str(output_path))
                mlflow.log_param("folder_read_pattern", folder_read_pattern)
                mlflow.log_param("file_read_pattern", file_read_pattern)
                mlflow.log_param("local_changes_hash", setuptools_scm.get_version())

                logging.info(
                    f"Starting MLflow experiment for dataset {dataset_name} session {session}..."
                )

                # Mock processing for different runs within the experiment
                for i in range(1, 11):  # 10 runs with varying parameters
                    # Start a child run under the main dataset-session run
                    with mlflow.start_run(nested=True):    

                        # Mock metric calculation                    
                        metric_measured = np.mean(data) * i 

                        # Log parameters and metrics specific to this run
                        mlflow.log_param("data_size", f"{i * 10}x100")
                        mlflow.log_param("run_iteration", i)
                        mlflow.log_param("run_id", mlflow.active_run().info.run_id)
                        mlflow.log_metric("metric_measured", metric_measured)

                        # Log the generated data as an artifact if desired
                        # Here, simulate an image or data file save path
                        image_path = writer.save_image(
                            image=data,
                            dataset_name=dataset_name,
                            session_number=session,
                            filename=f"image_run_{i}",
                        )
                        
                        mlflow.log_artifact(
                            image_path,
                            artifact_path=f"{dataset_name}/session_{session}/run_{i}",
                        )

                        logging.info(
                            f"Completed MLflow run iteration {i} for dataset {dataset_name} session {session}"
                        )

                logging.info(
                    f"Completed MLflow experiment for dataset {dataset_name} session {session}"
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
    parser.add_argument(
        "--experiment_name",
        type=str,
        help="Name of the experiment.",
        default="pipeline_test",
    )

    args = parser.parse_args()

    main(
        args.raw_data_path,
        args.output_path,
        args.folder_read_pattern,
        args.file_read_pattern,
        args.experiment_name,
    )
