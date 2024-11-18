import datetime
import logging
from pathlib import Path
from typing import List

import mlflow
import numpy as np
import setuptools_scm

from calcium_imaging_automation.core.reader import ReadAquiredData
from calcium_imaging_automation.core.writer import DatashuttleWrapper


def pipeline(
    raw_data_path: Path,
    output_path: Path,
    folder_read_pattern: str,
    file_read_pattern: List[str],
    experiment_name: str = "pipeline_test",
):
    # --- Setup logging and MLflow ---
    logging_setup(output_path)
    mlflow_setup(output_path)

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

    # --- Start processing ---
    for dataset in reader.datasets_paths:
        dataset_name = dataset.stem

        for session in range(0, number_of_tiffs):
            mlflow_set_experiment(experiment_name, dataset_name, session)

            # Generate mock data
            data = np.random.rand(100, 100)

            # Start a new MLflow experiment for each dataset-session
            with mlflow.start_run():  # this is the parent run
                mlflow_parent_run_logs(
                    dataset_name,
                    session,
                    raw_data_path,
                    output_path,
                    folder_read_pattern,
                    file_read_pattern,
                )

                logging.info(
                    f"Starting MLflow experiment for dataset {dataset_name} "
                    + f"session {session}..."
                )

                # Mock processing for different runs within the experiment
                for i in range(0, 10):  # n runs with varying parameters
                    # Start a child run under the main dataset-session run
                    with mlflow.start_run(nested=True):
                        # Mock metric calculation
                        metric_measured = np.mean(data) * i

                        # Log the generated data as an artifact if desired
                        # Here, simulate an image or data file save path
                        image_path = writer.save_image(
                            image=data,
                            dataset_name=dataset_name,
                            session_number=session,
                            filename=f"image_{mlflow.active_run().info.run_id}.png",
                        )

                        mlflow_log_run(
                            i,
                            dataset_name,
                            session,
                            metric_measured,
                            image_path,
                        )

                        logging.info(
                            f"Completed MLflow run iteration {i} for dataset "
                            + f"{dataset_name} session {session}"
                        )

                logging.info(
                    f"Completed MLflow experiment for dataset {dataset_name}"
                    + f" session {session}"
                )

    logging.info("Pipeline finished.")


def logging_setup(output_path: Path):
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


def mlflow_setup(output_path: Path):
    # --- Setup MLflow tracking ---
    mlflow_tracking_dir = output_path / "mlflow"
    mlflow.set_tracking_uri(str(mlflow_tracking_dir))


def mlflow_set_experiment(
    experiment_name: str, dataset_name: str, session: int
):
    # Start a new MLflow experiment for each dataset and session
    mlflow.set_experiment(
        f"{experiment_name}/{dataset_name}/session_{session}"
    )


def mlflow_parent_run_logs(
    dataset_name: str,
    session: int,
    raw_data_path: Path,
    output_path: Path,
    folder_read_pattern: str,
    file_read_pattern: List[str],
):
    # give specific name to the parent run
    mlflow.set_tag("mlflow.runName", f"{dataset_name}_session_{session}")

    # Log session-specific parameters
    mlflow.log_param("mlflow.Dataset", dataset_name)
    mlflow.log_param("session_number", session)
    mlflow.log_param("raw_data_path", str(raw_data_path))
    mlflow.log_param("output_path", str(output_path))
    mlflow.log_param("folder_read_pattern", folder_read_pattern)
    mlflow.log_param("file_read_pattern", file_read_pattern)
    mlflow.log_param("local_changes_hash", setuptools_scm.get_version())


def mlflow_log_run(
    i: int,
    dataset_name: str,
    session: int,
    metric_measured: float,
    image_path: Path,
):
    # give specific name to the run
    mlflow.set_tag("mlflow.runName", f"param_{i}")

    # Log parameters and metrics specific to this run
    mlflow.log_param("data_size", f"{i * 10}x100")
    mlflow.log_param("run_iteration", i)
    mlflow.log_param("run_id", mlflow.active_run().info.run_id)
    mlflow.log_metric("metric_measured", metric_measured)

    mlflow.log_artifact(
        # where I am storing the image according to Neuroblueprint
        # I think it gets copied in the mlflow data structure
        image_path,
        artifact_path=f"{dataset_name}/session_{session}/run_{i}",
    )
