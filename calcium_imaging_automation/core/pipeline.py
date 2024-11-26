import datetime
import logging
import time
from pathlib import Path
from typing import Callable, List

import mlflow
import setuptools_scm
import submitit
from submitit import AutoExecutor

from calcium_imaging_automation.core.reader import ReadAquiredData
from calcium_imaging_automation.core.writer import DatashuttleWrapper


def mlflow_orchestrator(
    raw_data_path: Path,
    output_path: Path,
    folder_read_pattern: str,
    file_read_pattern: List[str],
    preprocessing_function: Callable,
    compute_metric: Callable,
    experiment_name: str = "pipeline_test",
):
    # --- Setup logging and MLflow ---
    logging_setup(output_path)
    mlflow_setup(output_path)

    #  mkdir for submitit logs submitit / timestamp
    (output_path / "submitit").mkdir(exist_ok=True)

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
    results, errors = launch_job_array(
        datasets=reader.datasets_paths,
        output_path=output_path,
        analysis_pipeline=analysis_pipeline,
        writer=writer,
        preprocessing_function=preprocessing_function,
        compute_metric=compute_metric,
    )

    # --- Log all results with MLflow ---
    for dataset, result, error in zip(reader.dataset_names, results, errors):
        mlflow_set_experiment(experiment_name, dataset, 0)

        with mlflow.start_run():
            mlflow_parent_run_logs(
                dataset,
                0,
                raw_data_path,
                output_path,
                folder_read_pattern,
                file_read_pattern,
            )

            #  log error if any
            if error:
                mlflow.log_param("error", error)

            if result:
                mlflow.log_metric("stability", result)

            mlflow.end_run()

    logging.info("Pipeline finished.")


def launch_job_array(
    datasets,
    output_path,
    analysis_pipeline,
    writer,
    preprocessing_function,
    compute_metric,
):
    executor = AutoExecutor(folder=output_path / "submitit")
    executor.update_parameters(
        timeout_min=30,
        slurm_partition="fast",
        cpus_per_task=1,
        tasks_per_node=1,
        slurm_mem="16G",
        slurm_array_parallelism=20,
    )

    logging.info(f"Running {len(datasets)} jobs.")
    jobs = executor.map_array(
        analysis_pipeline,
        datasets,
        [writer.get_dataset_path(dataset.stem) for dataset in datasets],
        [preprocessing_function] * len(datasets),
        [compute_metric] * len(datasets),
    )

    results = []
    errors = []
    for job in jobs:
        while not job.done():
            time.sleep(10)
        try:
            results.append(job.result())
            errors.append(None)
        except submitit.core.utils.FailedJobError as e:
            logging.error(f"Job {job.job_id} failed: {e}")
            results.append(None)
            errors.append(job.stderr())

    return results, errors


def analysis_pipeline(
    dataset, output_path_dataset, preprocessing_function, compute_metric
):
    import os

    os.system("module load miniconda")
    os.system("source activate /nfs/nhome/live/lporta/.conda/envs/cimat")
    output_path_dataset = output_path_dataset / "ses-0/funcimg/"
    data = preprocessing_function(dataset, output_path_dataset)
    metric_measured = compute_metric(data)
    return metric_measured


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
    # image_path: Path,
):
    # give specific name to the run
    mlflow.set_tag("mlflow.runName", f"param_{i}")

    # Log parameters and metrics specific to this run
    mlflow.log_param("data_size", f"{i * 10}x100")
    mlflow.log_param("run_iteration", i)
    mlflow.log_param("run_id", mlflow.active_run().info.run_id)
    mlflow.log_metric("stability", metric_measured)

    # mlflow.log_artifact(
    #     # where I am storing the image according to Neuroblueprint
    #     # I think it gets copied in the mlflow data structure
    #     image_path,
    #     artifact_path=f"{dataset_name}/session_{session}/run_{i}",
    # )
