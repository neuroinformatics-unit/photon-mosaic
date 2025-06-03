import argparse
import importlib.resources as pkg_resources
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Run the photon-mosaic Snakemake pipeline."
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Path to your config.yaml file (defaults to internal one)",
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

    snakefile_path = pkg_resources.files("photon_mosaic").joinpath(
        "workflow", "Snakefile"
    )

    if args.config is None:
        config_path = pkg_resources.files("photon_mosaic").joinpath(
            "workflow", "config.yaml"
        )
    else:
        config_path = Path(args.config)

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

    subprocess.run(cmd)
