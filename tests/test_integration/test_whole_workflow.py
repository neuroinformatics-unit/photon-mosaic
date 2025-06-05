import shutil
import subprocess

import pytest
import yaml

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
                        "glob_naming_pattern_tif": ["2p_example_V1.tif"]
                    },
                }
            ],
            "output_patterns": ["2p_example_V1.tif"],
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
            "--printshellcmds",
            "--verbose",
            "--debug-dag",
        ],
        cwd=snake_test_env["workdir"],
        capture_output=True,
        text=True,
    )

    # Print detailed information about the failure
    print("\n=== Snakemake Dry Run Output ===")
    print("STDOUT:")
    print(result.stdout)
    print("\nSTDERR:")
    print(result.stderr)
    print("\nReturn Code:", result.returncode)
    print("=== End of Snakemake Dry Run Output ===\n")

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
            "--printshellcmds",
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

    # Print detailed information about the execution
    print("\n=== Snakemake Execution Output ===")
    print("STDOUT:")
    print(result.stdout)
    print("\nSTDERR:")
    print(result.stderr)
    print("\nReturn Code:", result.returncode)
    print("=== End of Snakemake Execution Output ===\n")

    assert result.returncode == 0, (
        f"Snakemake execution failed:\nSTDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    output_base = (
        snake_test_env["workdir"]
        / "processed"
        / "sub-0_001"
        / "ses-0"
        / "funcimg"
        / "suite2p"
        / "plane0"
    )

    # Print information about the expected output files
    print("\n=== Expected Output Files ===")
    print(f"Checking for files in: {output_base}")
    print(
        "Directory contents:",
        list(output_base.iterdir())
        if output_base.exists()
        else "Directory does not exist",
    )
    print("=== End of Expected Output Files ===\n")

    assert (output_base / "F.npy").exists(), "Missing output: F.npy"
    assert (output_base / "data.bin").exists(), "Missing output: data.bin"


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
