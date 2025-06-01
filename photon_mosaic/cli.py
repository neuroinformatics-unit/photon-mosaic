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
        "--configfile",
        str(config_path),
        "--jobs",
        args.jobs,
    ]

    if args.dry_run:
        cmd.append("--dry-run")

    if args.extra:
        cmd.extend(args.extra)

    subprocess.run(cmd)
