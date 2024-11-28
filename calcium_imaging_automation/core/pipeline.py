import datetime
import logging
import time
from pathlib import Path
from typing import Callable, List

import pandas as pd
import submitit
from submitit import AutoExecutor

from calcium_imaging_automation.core.reader import ReadAquiredData
from calcium_imaging_automation.core.writer import DatashuttleWrapper


def orchestrator(
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

    # save the results and errors as csv
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path / "results.csv")
    errors_df = pd.DataFrame(errors)
    errors_df.to_csv(output_path / "errors.csv")

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
    try:
        data = preprocessing_function(dataset, output_path_dataset)
        metric_measured = compute_metric(data)
        with open(output_path_dataset / "metric.txt", "w") as f:
            f.write(str(metric_measured))
    except Exception as e:
        with open(output_path_dataset / "error.txt", "w") as f:
            f.write(str(e.args))
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
