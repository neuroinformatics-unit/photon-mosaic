import subprocess

import numpy as np
import yaml
from tifffile import imread
from pathlib import Path

from photon_mosaic import get_snakefile_path


def run_snakemake(workdir, configfile, dry_run=False):
    """Helper function to run snakemake with common parameters."""
    snakefile = str(get_snakefile_path())
    cmd = [
        "snakemake",
        "--cores",
        "1",
        "--verbose",
        "--keep-going",
        "-s",
        str(snakefile),
        "--configfile",
        str(configfile),
        "--debug-dag",
    ]

    if dry_run:
        cmd.insert(1, "--dry-run")

    result = subprocess.run(
        cmd,
        cwd=workdir,
        capture_output=True,
        text=True,
    )

    return result


def check_output_files(workdir, datasets, tiffs, check_enhanced=False):
    """Helper function to check output files."""
    for sub_idx, dataset in enumerate(datasets):
        for ses_idx, tiff in enumerate(tiffs):
            output_base = (
                workdir
                / "processed"
                / f"sub-{sub_idx}_{dataset}"
                / f"ses-{ses_idx}"
                / "funcimg"
                / "suite2p"
                / "plane0"
            )

            # Print information about the expected output files
            print(
                f"\n=== Checking files for {dataset}/ses-{ses_idx}/{tiff} ==="
            )
            print(f"Checking for files in: {output_base}")
            print(
                "Directory contents:",
                list(output_base.iterdir())
                if output_base.exists()
                else "Directory does not exist",
            )
            print("=== End of Expected Output Files ===\n")

            assert (
                output_base / "F.npy"
            ).exists(), (
                f"Missing output: F.npy for {dataset}/ses-{ses_idx}/{tiff}"
            )
            assert (
                output_base / "data.bin"
            ).exists(), (
                f"Missing output: data.bin for {dataset}/ses-{ses_idx}/{tiff}"
            )

            if check_enhanced:
                enhanced_file = (
                    workdir
                    / "processed"
                    / f"sub-{sub_idx}_{dataset}"
                    / f"ses-{ses_idx}"
                    / "funcimg"
                    / f"enhanced_2p_example_V1_{ses_idx+1:02d}.tif"
                )

                assert (
                    enhanced_file.exists()
                ), f"Missing enhanced output: {enhanced_file}"

                # Verify that the enhanced file contains valid image data
                enhanced_data = imread(enhanced_file)
                assert (
                    enhanced_data.dtype == np.int16
                ), f"Enhanced file {enhanced_file} has "
                f"incorrect data type: {enhanced_data.dtype}"
                assert (
                    enhanced_data.size > 0
                ), f"Enhanced file {enhanced_file} is empty"


def test_snakemake_dry_run(snake_test_env):
    """Test that snakemake can do a dry run."""
    result = run_snakemake(
        snake_test_env["workdir"], snake_test_env["configfile"], dry_run=True
    )
    assert result.returncode == 0, (
        f"Snakemake dry-run failed:\nSTDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )


def test_snakemake_execution(snake_test_env):
    """Test that snakemake can execute the workflow."""
    result = run_snakemake(
        snake_test_env["workdir"], snake_test_env["configfile"]
    )
    assert result.returncode == 0, (
        f"Snakemake execution failed:\nSTDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    # Check that output files were created for each dataset and tiff
    datasets = ["001", "002", "003"]
    tiffs = ["2p_example_V1_01.tif", "2p_example_V1_02.tif"]
    check_output_files(snake_test_env["workdir"], datasets, tiffs)


def test_snakemake_with_contrast(snake_test_env, test_config_with_contrast):
    """
    Test that snakemake can execute the workflow with contrast enhancement
    preprocessing.
    """
    # Update config with contrast enhancement settings
    config = test_config_with_contrast.copy()
    config["raw_data_base"] = str(Path(snake_test_env["workdir"]) / "raw")
    config["processed_data_base"] = str(
        Path(snake_test_env["workdir"]) / "processed"
    )

    # Create config file
    config_path = Path(snake_test_env["workdir"]) / "config.yaml"
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, default_style='"', allow_unicode=True)

    # Run snakemake with real contrast enhancement
    result = run_snakemake(snake_test_env["workdir"], config_path)
    assert result.returncode == 0, (
        f"Snakemake execution with contrast enhancement failed:\nSTDOUT:\n"
        f"{result.stdout}\nSTDERR:\n{result.stderr}"
    )

    # Check that enhanced output files were created and contain valid image
    # data
    datasets = ["001", "002", "003"]
    tiffs = ["2p_example_V1_01.tif", "2p_example_V1_02.tif"]
    check_output_files(
        snake_test_env["workdir"], datasets, tiffs, check_enhanced=True
    )

 
def test_photon_mosaic_cli_dry_run(snake_test_env):
    """Test that photon-mosaic can do a dry run."""
    result = run_snakemake(
        snake_test_env["workdir"], snake_test_env["configfile"], dry_run=True
    )
    assert result.returncode == 0, (
        f"photon-mosaic CLI dry-run failed:\nSTDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )
