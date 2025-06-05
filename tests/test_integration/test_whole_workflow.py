import shutil
import subprocess

import numpy as np
import pytest
import yaml
from tifffile import imread
from pathlib import Path

from photon_mosaic import get_snakefile_path


@pytest.fixture
def test_config():
    """Create a test configuration."""
    return {
        "raw_data_base": "raw",  # This will be replaced in snake_test_env
        "processed_data_base": "processed",  # This will be replaced in
        # snake_test_env
        "dataset_discovery": {
            "pattern": "^.*$",  # Match all directories for testing
            "exclude_patterns": [],  # Don't exclude anything
        },
        "preprocessing": {
            "steps": [
                {
                    "name": "noop",
                    "kwargs": {
                        "glob_naming_pattern_tif": [
                            "2p_example_V1_01.tif",
                            "2p_example_V1_02.tif",
                        ]
                    },
                }
            ],
            "output_patterns": [
                "2p_example_V1_01.tif",
                "2p_example_V1_02.tif",
            ],
        },
        "suite2p_ops": {
            "fs": 6.0,
            "nplanes": 1,
            "tau": 0.7,
            "nonrigid": True,
            "diameter": 8,
        },
        "use_slurm": False,
    }


@pytest.fixture
def test_config_with_contrast():
    """
    Create a test configuration with contrast enhancement preprocessing step.
    """
    return {
        "raw_data_base": "raw",  # This will be replaced in snake_test_env
        "processed_data_base": "processed",  # This will be replaced in
        # snake_test_env
        "dataset_discovery": {
            "pattern": "^.*$",  # Match all directories for testing
            "exclude_patterns": [],  # Don't exclude anything
        },
        "preprocessing": {
            "steps": [
                {
                    "name": "contrast",
                    "kwargs": {
                        "clip_limit": 2.0,
                        "glob_naming_pattern_tif": [
                            "2p_example_V1_01.tif",
                            "2p_example_V1_02.tif",
                        ],
                    },
                }
            ],
            "output_patterns": [
                "enhanced_2p_example_V1_01.tif",
                "enhanced_2p_example_V1_02.tif",
            ],
        },
        "suite2p_ops": {
            "fs": 6.0,
            "nplanes": 1,
            "tau": 0.7,
            "nonrigid": True,
            "diameter": 8,
        },
        "use_slurm": False,
    }


@pytest.fixture
def snake_test_env(tmp_path, test_config):
    """
    Fixture that sets up the test environment with data and configuration.
    """
    raw_data = tmp_path / "raw"
    shutil.copytree("tests/data", raw_data)
    processed_data = tmp_path / "processed"
    processed_data.mkdir()

    # Update paths in config
    config = test_config.copy()
    config["raw_data_base"] = str(raw_data)
    config["processed_data_base"] = str(processed_data)

    # Create config file
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, default_style='"', allow_unicode=True)

    return {
        "workdir": tmp_path,
        "configfile": config_path,
    }


def test_snakemake_dry_run(snake_test_env):
    """Test that snakemake can do a dry run."""
    snakefile = str(get_snakefile_path())

    result = subprocess.run(
        [
            "snakemake",
            "--dry-run",
            "-s",
            snakefile,
            "--configfile",
            str(snake_test_env["configfile"]),
            "--verbose",
            "--debug-dag",
        ],
        cwd=snake_test_env["workdir"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"Snakemake dry-run failed:\nSTDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )


def test_snakemake_execution(snake_test_env):
    """Test that snakemake can execute the workflow."""
    snakefile = str(get_snakefile_path())

    result = subprocess.run(
        [
            "snakemake",
            "--cores",
            "1",
            "--verbose",
            "--keep-going",
            "-s",
            snakefile,
            "--configfile",
            str(snake_test_env["configfile"]),
            "--debug-dag",  # Add DAG debugging
        ],
        cwd=snake_test_env["workdir"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"Snakemake execution failed:\nSTDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    # Check that output files were created for each dataset and tiff
    datasets = ["001", "002", "003"]
    tiffs = ["2p_example_V1_01.tif", "2p_example_V1_02.tif"]

    for sub_idx, dataset in enumerate(datasets):
        for ses_idx, tiff in enumerate(tiffs):
            output_base = (
                snake_test_env["workdir"]
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


def test_snakemake_with_contrast(snake_test_env, test_config_with_contrast):
    """Test that snakemake can execute the workflow with contrast enhancement
    preprocessing."""
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

    snakefile = str(get_snakefile_path())

    # Run snakemake with real contrast enhancement
    result = subprocess.run(
        [
            "snakemake",
            "--cores",
            "1",
            "--verbose",
            "--keep-going",
            "-s",
            str(snakefile),
            "--configfile",
            str(config_path),
            "--debug-dag",
        ],
        cwd=snake_test_env["workdir"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"Snakemake execution with contrast enhancement failed:\nSTDOUT:\n"
        f"{result.stdout}\nSTDERR:\n{result.stderr}"
    )

    # Check that enhanced output files were created and contain valid image
    # data
    datasets = ["001", "002", "003"]
    for sub_idx, dataset in enumerate(datasets):
        for ses_idx in range(2):
            output_base = (
                snake_test_env["workdir"]
                / "processed"
                / f"sub-{sub_idx}_{dataset}"
                / f"ses-{ses_idx}"
                / "funcimg"
                / "suite2p"
                / "plane0"
            )

            # Check for enhanced files
            enhanced_file = (
                snake_test_env["workdir"]
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

            assert (
                output_base / "F.npy"
            ).exists(), f"Missing output: F.npy for {dataset}/ses-{ses_idx}"
            assert (
                output_base / "data.bin"
            ).exists(), f"Missing output: data.bin for {dataset}/ses-{ses_idx}"


def test_photon_mosaic_cli_dry_run(snake_test_env):
    result = subprocess.run(
        [
            "photon-mosaic",
            "--config",
            str(snake_test_env["configfile"]),
            "--jobs",
            "1",
            "--dry-run",
        ],
        cwd=snake_test_env["workdir"],
        capture_output=True,
        text=True,
    )

    print(result.stdout)
    assert (
        result.returncode == 0
    ), f"photon-mosaic CLI dry-run failed:\n{result.stderr}"
