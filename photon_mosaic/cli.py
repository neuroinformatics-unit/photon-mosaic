import argparse
import importlib.resources as pkg_resources
import subprocess
from datetime import datetime
from pathlib import Path

import yaml


def main():
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

    # Ensure ~/.photon_mosaic/config.yaml exists
    user_config_dir = Path.home() / ".photon_mosaic"
    user_config_path = user_config_dir / "config.yaml"
    if not user_config_path.exists():
        user_config_dir.mkdir(parents=True, exist_ok=True)
        default_config_path = pkg_resources.files("photon_mosaic").joinpath(
            "workflow", "config.yaml"
        )
        with (
            open(default_config_path, "r") as src,
            open(user_config_path, "w") as dst,
        ):
            dst.write(src.read())

    # Determine which config to use
    if args.config is not None:
        config_path = Path(args.config)
    elif user_config_path.exists():
        config_path = user_config_path
    else:
        config_path = pkg_resources.files("photon_mosaic").joinpath(
            "workflow", "config.yaml"
        )

    # Load config
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Apply CLI overrides
    if args.raw_data_base is not None:
        config["raw_data_base"] = args.raw_data_base
    if args.processed_data_base is not None:
        config["processed_data_base"] = args.processed_data_base

    # Create photon-mosaic directory with logs and configs subdirectories
    photon_mosaic_dir = Path("photon-mosaic")
    logs_dir = photon_mosaic_dir / "logs"
    configs_dir = photon_mosaic_dir / "configs"
    photon_mosaic_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    configs_dir.mkdir(exist_ok=True)

    # Generate timestamp for this run
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Save config with timestamp
    config_filename = f"config_{timestamp}.yaml"
    config_path = configs_dir / config_filename
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    snakefile_path = pkg_resources.files("photon_mosaic").joinpath(
        "workflow", "Snakefile"
    )

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
