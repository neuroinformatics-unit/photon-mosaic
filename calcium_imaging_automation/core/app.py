import argparse
from pathlib import Path

from calcium_imaging_automation.core.pipeline import pipeline

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

    pipeline(
        args.raw_data_path,
        args.output_path,
        args.folder_read_pattern,
        args.file_read_pattern,
        args.experiment_name,
    )
