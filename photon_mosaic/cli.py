"""
Command line interface for photon-mosaic.
"""

import argparse
import importlib.resources as pkg_resources
import logging
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

    Command Line Arguments
    ---------------------
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
    --log-level : str, default="INFO"
        Log level.
    --verbose, -v : flag, optional
        Enable verbose output (default is quiet mode).
    --reset-config : flag, optional
        Reset the config file to the default values.
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
        default=None,
        help="Override raw_data_base in config file",
    )
    parser.add_argument(
        "--processed_data_base",
        default=None,
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
        "--log-level",
        default="INFO",
        help="Log level",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose output (default is quiet mode)",
    )
    parser.add_argument(
        "extra",
        nargs=argparse.REMAINDER,
        help="Additional arguments to snakemake",
    )
    parser.add_argument(
        "--reset-config",
        action="store_true",
        help="Reset the config file to the default values",
    )

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level)
    logger = logging.getLogger(__name__)
    logger.debug("Starting photon-mosaic CLI")

    # Ensure ~/.photon_mosaic/config.yaml exists, if not, create it
    default_config_dif = Path.home() / ".photon_mosaic"
    default_config_path = default_config_dif / "config.yaml"
    if not default_config_path.exists() or args.reset_config:
        logger.debug("Creating default config file")
        default_config_dif.mkdir(parents=True, exist_ok=True)
        source_config_path = pkg_resources.files("photon_mosaic").joinpath(
            "workflow", "config.yaml"
        )
        with (
            source_config_path.open("rb") as src,
            open(default_config_path, "wb") as dst,
        ):
            dst.write(src.read())

    # Determine which config to use
    if args.config is not None:
        logger.debug(f"Using config file: {args.config}")
        # Take the path provided by the user
        config_path = Path(args.config)
    elif default_config_path.exists():
        logger.debug("Using default config file")
        # Use the default config
        config_path = default_config_path

    # Load config
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Apply CLI overrides
    if args.raw_data_base is not None:
        logger.debug(f"Overriding raw_data_base to: {args.raw_data_base}")
        config["raw_data_base"] = args.raw_data_base
    else:
        logger.debug(
            f"Using raw_data_base from config file: {config['raw_data_base']}"
        )

    if args.processed_data_base is not None:
        logger.debug(
            f"Overriding processed_data_base to: {args.processed_data_base}"
        )
        config["processed_data_base"] = args.processed_data_base
    else:
        logger.debug(
            "Using processed_data_base from config file: "
            f"{config['processed_data_base']}"
        )

    raw_data_base = Path(config["raw_data_base"]).resolve()

    # Append derivatives to the processed_data_base if it doesn't end with
    # /derivatives
    processed_data_base = Path(config["processed_data_base"]).resolve()
    if processed_data_base.name != "derivatives":
        processed_data_base = processed_data_base / "derivatives"
    config["processed_data_base"] = str(processed_data_base)

    # Only update the default config file if it's the one being used
    # Don't modify the user's default config when
    # they've specified a custom config
    if args.config is None:
        # Change the values of processed_data_base in the config file saved in
        # the .photon_mosaic/config.yaml without changing the other values and
        # losing the comments
        with open(default_config_path, "r") as f:
            config_lines = f.readlines()
        with open(default_config_path, "w") as f:
            for line in config_lines:
                if line.startswith("processed_data_base:"):
                    f.write(f"processed_data_base: {processed_data_base}\n")
                elif line.startswith("raw_data_base:"):
                    f.write(f"raw_data_base: {raw_data_base}\n")
                else:
                    f.write(line)

    # Create photon-mosaic directory with logs and configs subdirectories
    output_dir = processed_data_base / "photon-mosaic"
    logger.debug(f"Creating output directory: {output_dir}")
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
        logger.debug(f"Saving config to: {config_path}")
        yaml.dump(config, f)

    snakefile_path = get_snakefile_path()

    logger.debug(f"Launching snakemake with snakefile: {snakefile_path}")
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
    if not args.verbose:
        # Default to quiet mode unless --verbose is specified
        cmd.append("--quiet")
    if config["use_slurm"] == "slurm":
        cmd.extend(["--executor", "slurm"])
    if args.extra:
        cmd.extend(args.extra)

    # Save logs with timestamp
    log_filename = f"snakemake_{timestamp}.log"
    log_path = logs_dir / log_filename
    with open(log_path, "w") as logfile:
        logger.debug(f"Saving logs to: {log_path}")
        logger.info(f"Launching snakemake with command: {' '.join(cmd)}")
        result = subprocess.run(cmd, stdout=logfile, stderr=logfile)

        # If workflow is locked, automatically unlock and retry
        if result.returncode != 0 and not args.unlock:
            # Check if the error is related to locking by reading the log
            with open(log_path, "r") as log_read:
                log_content = log_read.read()
                if (
                    "locked" in log_content.lower()
                    or "lock" in log_content.lower()
                ):
                    logger.info(
                        "Directory appears to be locked. "
                        "Unlocking automatically and retrying..."
                    )

                    # Create unlock command
                    unlock_cmd = cmd.copy()
                    unlock_cmd.append("--unlock")

                    # Run unlock command
                    with open(log_path, "a") as logfile_append:
                        logfile_append.write("\n--- Auto-unlock attempt ---\n")
                        unlock_result = subprocess.run(
                            unlock_cmd,
                            stdout=logfile_append,
                            stderr=logfile_append,
                        )

                    # If unlock succeeded, retry the original command
                    if unlock_result.returncode == 0:
                        logger.info(
                            "Unlock successful. Retrying original command..."
                        )
                        with open(log_path, "a") as logfile_append:
                            logfile_append.write(
                                "\n--- Retry after unlock ---\n"
                            )
                            result = subprocess.run(
                                cmd,
                                stdout=logfile_append,
                                stderr=logfile_append,
                            )

        if result.returncode == 0:
            logging.info("Snakemake pipeline completed successfully.")
        else:
            logging.info(
                "Snakemake pipeline failed with exit code "
                f"{result.returncode}. Check the log file at "
                f"{log_path} for details."
            )
