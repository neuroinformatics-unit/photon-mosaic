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


def create_argument_parser():
    """Create and configure the argument parser for the CLI.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser with all CLI options
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
    return parser


def ensure_default_config(reset_config=False):
    """Ensure the default config file exists, create if missing or
    reset requested.

    Parameters
    ----------
    reset_config : bool
        Whether to reset the config file to defaults

    Returns
    -------
    Path
        Path to the default config file
    """
    logger = logging.getLogger(__name__)
    default_config_dir = Path.home() / ".photon_mosaic"
    default_config_path = default_config_dir / "config.yaml"

    if not default_config_path.exists() or reset_config:
        logger.debug("Creating default config file")
        default_config_dir.mkdir(parents=True, exist_ok=True)
        source_config_path = pkg_resources.files("photon_mosaic").joinpath(
            "workflow", "config.yaml"
        )
        with (
            source_config_path.open("rb") as src,
            open(default_config_path, "wb") as dst,
        ):
            dst.write(src.read())

    return default_config_path


def load_and_process_config(args):
    """Load configuration file and apply CLI overrides.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line arguments

    Returns
    -------
    tuple[dict, Path]
        Configuration dictionary and path to config file used
    """
    logger = logging.getLogger(__name__)

    # Ensure default config exists
    default_config_path = ensure_default_config(args.reset_config)

    # Determine which config to use
    if args.config is not None:
        logger.debug(f"Using config file: {args.config}")
        config_path = Path(args.config)
    else:
        logger.debug("Using default config file")
        config_path = default_config_path

    # Load config
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    logger.debug(f"Loaded config from {config_path}:")
    logger.debug(f"use_slurm: {config.get('use_slurm', 'NOT SET')}")
    logger.debug(f"slurm config: {config.get('slurm', 'NOT SET')}")

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

    # Process paths
    raw_data_base = Path(config["raw_data_base"]).resolve()
    processed_data_base = Path(config["processed_data_base"]).resolve()

    # Append derivatives to the processed_data_base if
    # it doesn't end with /derivatives
    if processed_data_base.name != "derivatives":
        processed_data_base = processed_data_base / "derivatives"
    config["processed_data_base"] = str(processed_data_base)

    # Update default config file if it's the one being used
    if args.config is None:
        update_default_config(
            default_config_path, raw_data_base, processed_data_base
        )

    return config, config_path


def update_default_config(config_path, raw_data_base, processed_data_base):
    """Update the default config file with new base paths while
    preserving comments.

    Parameters
    ----------
    config_path : Path
        Path to the config file to update
    raw_data_base : Path
        Raw data base path
    processed_data_base : Path
        Processed data base path
    """
    with open(config_path, "r") as f:
        config_lines = f.readlines()
    with open(config_path, "w") as f:
        for line in config_lines:
            if line.startswith("processed_data_base:"):
                f.write(f"processed_data_base: {processed_data_base}\n")
            elif line.startswith("raw_data_base:"):
                f.write(f"raw_data_base: {raw_data_base}\n")
            else:
                f.write(line)


def setup_output_directories(processed_data_base):
    """Create necessary output directories for the pipeline.

    Parameters
    ----------
    processed_data_base : Path
        Base path for processed data

    Returns
    -------
    tuple[Path, Path, Path]
        Paths to output directory, logs directory, and configs directory
    """
    logger = logging.getLogger(__name__)

    output_dir = processed_data_base / "photon-mosaic"
    logger.debug(f"Creating output directory: {output_dir}")

    logs_dir = output_dir / "logs"
    configs_dir = output_dir / "configs"

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    configs_dir.mkdir(exist_ok=True)

    return output_dir, logs_dir, configs_dir


def save_timestamped_config(config, configs_dir):
    """Save configuration with timestamp for reproducibility.

    Parameters
    ----------
    config : dict
        Configuration dictionary to save
    configs_dir : Path
        Directory to save config files

    Returns
    -------
    tuple[str, Path]
        Timestamp string and path to saved config file
    """
    logger = logging.getLogger(__name__)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    config_filename = f"config_{timestamp}.yaml"
    config_path = configs_dir / config_filename

    with open(config_path, "w") as f:
        logger.debug(f"Saving config to: {config_path}")
        yaml.dump(config, f)

    return timestamp, config_path


def build_snakemake_command(args, config_path):
    """Build the base snakemake command with common arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line arguments
    config_path : Path
        Path to the config file to use

    Returns
    -------
    list[str]
        Snakemake command as list of strings
    """
    logger = logging.getLogger(__name__)

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
    if args.latency_wait:
        cmd.extend(["--latency-wait", str(args.latency_wait)])
    if not args.verbose:
        # Default to quiet mode unless --verbose is specified
        cmd.append("--quiet")

    return cmd


def configure_slurm_execution(cmd, config):
    """Configure SLURM execution options if enabled.

    Parameters
    ----------
    cmd : list[str]
        Base snakemake command to extend
    config : dict
        Configuration dictionary

    Returns
    -------
    list[str]
        Updated command with SLURM configuration
    """
    logger = logging.getLogger(__name__)

    if not config.get("use_slurm", False):
        logger.info("SLURM execution disabled - running locally")
        return cmd

    logger.info("SLURM execution enabled - configuring SLURM executor")
    cmd.extend(["--executor", "slurm"])

    # Add SLURM-specific arguments
    slurm_config = config.get("slurm", {})
    logger.info(f"SLURM configuration loaded: {slurm_config}")

    # Resources that should NOT be passed via --default-resources
    # because they're already set at rule level and may cause conflicts
    # (particularly gpu-related resources that can trigger TRES conflicts)
    # Note: tasks_per_gpu is NOT excluded as it's needed
    # to control --ntasks-per-gpu behavior
    exclude_from_defaults = {"gpu", "gres", "cpus_per_gpu"}

    # Create default resources string for SLURM by iterating all provided keys
    default_resources = []
    for key, value in slurm_config.items():
        # Skip empty/null values
        if value is None:
            continue
        # Skip resources that should only be set at rule level
        if key in exclude_from_defaults:
            logger.info(
                f"Skipping {key} in --default-resources "
                f"(set at rule level to avoid conflicts): {value}"
            )
            continue
        # Quote string values so snakemake parses them correctly
        if isinstance(value, str):
            default_resources.append(f"{key}='{value}'")
        # Use lower-case for booleans (true/false)
        elif isinstance(value, bool):
            default_resources.append(f"{key}={str(value).lower()}")
        else:
            default_resources.append(f"{key}={value}")

        logger.info(f"Using SLURM {key}: {value}")

    if default_resources:
        # Join all default resources with commas and pass as single argument
        resources_str = ",".join(default_resources)
        cmd.extend(["--default-resources", resources_str])
        logger.info(f"SLURM default resources: {resources_str}")

    logger.info("SLURM executor configured successfully")
    return cmd


def execute_pipeline_with_retry(cmd, log_path):
    """Execute the snakemake pipeline with automatic
    unlock and retry on lock errors.

    Parameters
    ----------
    cmd : list[str]
        Snakemake command to execute
    log_path : Path
        Path to save execution logs

    Returns
    -------
    int
        Return code from the snakemake execution
    """
    logger = logging.getLogger(__name__)

    with open(log_path, "w") as logfile:
        logger.debug(f"Saving logs to: {log_path}")
        logger.info(f"Launching snakemake with command: {' '.join(cmd)}")
        result = subprocess.run(cmd, stdout=logfile, stderr=logfile)

        # If workflow is locked, automatically unlock and retry
        if result.returncode != 0:
            # Check if the error is related to locking by reading the log
            with open(log_path, "r") as log_read:
                log_content = log_read.read()
                if (
                    "locked" in log_content.lower()
                    or "lock" in log_content.lower()
                ):
                    logger.warning(
                        "Workflow appears to be locked. "
                        "Attempting automatic unlock and retry..."
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
                            "Automatic unlock successful. "
                            "Retrying pipeline execution..."
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
                    else:
                        logger.error(
                            "Automatic unlock failed. "
                            "Please check the log file for details."
                        )

    return result.returncode


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
    # Parse command line arguments
    parser = create_argument_parser()
    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(level=args.log_level)
    logger = logging.getLogger(__name__)
    logger.debug("Starting photon-mosaic CLI")

    # Load and process configuration
    config, _ = load_and_process_config(args)

    # Set up output directories
    processed_data_base = Path(config["processed_data_base"])
    output_dir, logs_dir, configs_dir = setup_output_directories(
        processed_data_base
    )

    # Save timestamped config for reproducibility
    timestamp, config_path = save_timestamped_config(config, configs_dir)

    # Build snakemake command
    cmd = build_snakemake_command(args, config_path)

    # Configure SLURM execution if enabled
    cmd = configure_slurm_execution(cmd, config)

    # Add extra arguments if provided
    if args.extra:
        cmd.extend(args.extra)

    # Execute pipeline with automatic retry on lock errors
    log_filename = f"snakemake_{timestamp}.log"
    log_path = logs_dir / log_filename
    return_code = execute_pipeline_with_retry(cmd, log_path)

    # Report final status
    if return_code == 0:
        logging.info("Snakemake pipeline completed successfully.")
    else:
        logging.info(
            f"Snakemake pipeline failed with exit code {return_code}. "
            f"Check the log file at {log_path} for details."
        )
