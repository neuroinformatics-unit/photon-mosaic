import argparse
import importlib.resources as pkg_resources
import subprocess
from datetime import datetime
from pathlib import Path

import yaml

from photon_mosaic import get_snakefile_path


def main():
    """Run the photon-mosaic Snakemake pipeline for automated and reproducible
     analysis of multiphoton calcium imaging datasets.

    This pipeline integrates Suite2p for image registration and signal
    extraction, with a standardized output folder structure following
    the NeuroBlueprint specification. It is designed for labs that
    store their data on servers connected to HPC clusters and want to
    batch-process multiple imaging sessions in parallel.

    Parameters
    ----------
    --config : str, optional
        Path to your config.yaml file. If not provided, uses
        ~/.photon_mosaic/config.yaml.
    --raw_data_base : str, optional
        Override raw_data_base in config file (path to raw imaging data).
    --processed_data_base : str, optional
        Override processed_data_base in config file (path for processed
        outputs).
    --jobs : str, default="1"
        Number of parallel jobs to run.
    --dry-run : bool, optional
        Perform a dry run to preview the workflow without executing.
    --forcerun : str, optional
        Force re-run of a specific rule (e.g., 'suite2p').
    --rerun-incomplete : bool, optional
        Rerun any incomplete jobs.
    --unlock : bool, optional
        Unlock the workflow if it's in a locked state.
    --latency-wait : int, default=10
        Time to wait before checking if output files are ready.
    --executor : str, default="slurm"
        Executor to use.
    extra : list
        Additional arguments to pass to snakemake.

    Notes
    -----
    The pipeline will:
    1. Create a timestamped config file in derivatives/photon-mosaic/configs/
    2. Save execution logs in derivatives/photon-mosaic/logs/
    3. Process all TIFF files found in the raw data directory
    4. Generate standardized outputs following NeuroBlueprint specification
    """
    parser = argparse.ArgumentParser(
        description="Run the photon-mosaic Snakemake pipeline."
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Path to your config.yaml file.",
    )
    parser.add_argument(
        "--raw_data_base",
        default=".",
        help="Override raw_data_base in config file",
    )
    parser.add_argument(
        "--processed_data_base",
        default=".",
        help="Override processed_data_base in config file",
    )
    parser.add_argument(
        "--jobs", default="1", help="Number of parallel jobs to run"
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Perform a dry run"
    )
    parser.add_argument(
        "--forcerun",
        default=None,
        help="Force re-run of a specific rule",
    )
    parser.add_argument(
        "--rerun-incomplete",
        action="store_true",
        help="Rerun any incomplete jobs",
    )
    parser.add_argument(
        "--unlock",
        action="store_true",
        help="Unlock the workflow if it's in a locked state",
    )
    parser.add_argument(
        "--latency-wait",
        default=10,
        help="Time to wait before checking if output files are ready",
    )
    parser.add_argument(
        "--executor",
        default="slurm",
        help="Executor to use",
    )
    parser.add_argument(
        "extra",
        nargs=argparse.REMAINDER,
        help="Additional arguments to snakemake",
    )

    args = parser.parse_args()

    # Ensure ~/.photon_mosaic/config.yaml exists, if not, create it
    default_config_dif = Path.home() / ".photon_mosaic"
    default_config_path = default_config_dif / "config.yaml"
    if not default_config_path.exists():
        default_config_dif.mkdir(parents=True, exist_ok=True)
        default_config_path = pkg_resources.files("photon_mosaic").joinpath(
            "workflow", "config.yaml"
        )
        with (
            open(default_config_path, "r") as src,
            open(default_config_path, "w") as dst,
        ):
            dst.write(src.read())

    # Determine which config to use
    if args.config is not None:
        # Take the path provided by the user
        config_path = Path(args.config)
    elif default_config_path.exists():
        # Use the default config
        config_path = default_config_path

    # Load config
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Apply CLI overrides
    if args.raw_data_base is not None:
        config["raw_data_base"] = args.raw_data_base
    if args.processed_data_base is not None:
        config["processed_data_base"] = args.processed_data_base

    # Create photon-mosaic directory with logs and configs subdirectories
    output_dir = (
        Path(args.processed_data_base) / "derivatives" / "photon-mosaic"
    )
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    logs_dir = output_dir / "logs"
    configs_dir = output_dir / "configs"
    output_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    configs_dir.mkdir(exist_ok=True)

    # Generate timestamp for this run
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Save config with timestamp
    config_filename = f"config_{timestamp}.yaml"
    config_path = configs_dir / config_filename
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    snakefile_path = get_snakefile_path()

    cmd = [
        "snakemake",
        "--snakefile",
        str(snakefile_path),
        "--jobs",
        str(args.jobs),
        "--configfile",
        str(config_path),
    ]

    if args.dry_run:
        cmd.append("--dry-run")
    if args.forcerun:
        cmd.extend(["--forcerun", args.forcerun])
    if args.rerun_incomplete:
        cmd.append("--rerun-incomplete")
    if args.unlock:
        cmd.append("--unlock")
    if args.latency_wait:
        cmd.extend(["--latency-wait", str(args.latency_wait)])
    if args.executor:
        cmd.extend(["--executor", args.executor])
    if args.extra:
        cmd.extend(args.extra)

    # Save logs with timestamp
    log_filename = f"snakemake_{timestamp}.log"
    log_path = logs_dir / log_filename
    with open(log_path, "w") as logfile:
        subprocess.run(cmd, stdout=logfile, stderr=logfile)
