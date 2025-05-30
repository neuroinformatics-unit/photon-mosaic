import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

WORKFLOW_PATH = Path(__file__).parent.parent.parent / "workflow" / "Snakefile"


@pytest.fixture
def snake_test_env(tmp_path):
    raw_data = tmp_path / "raw"
    shutil.copytree("tests/data", raw_data)

    processed_data = tmp_path / "processed"
    processed_data.mkdir()

    config = {
        "raw_data_base": str(raw_data),
        "processed_data_base": str(processed_data),
        "suite2p_ops": {
            "fs": 6.0,
            "nplanes": 1,
            "tau": 0.7,
            "nonrigid": True,
            "diameter": 8,
        },
        "slurm": {
            "use_slurm": False,
        },
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    return {
        "workdir": tmp_path,
        "configfile": config_path,
    }


def test_snakemake_dry_run(snake_test_env):
    result = subprocess.run(
        [
            "snakemake",
            "--dry-run",
            "-s",
            str(WORKFLOW_PATH),
            "--configfile",
            str(snake_test_env["configfile"]),
        ],
        cwd=snake_test_env["workdir"],
        capture_output=True,
        text=True,
    )

    print(result.stdout)
    assert (
        result.returncode == 0
    ), f"Snakemake dry-run failed:\n{result.stderr}"


def test_snakemake_execution(snake_test_env):
    result = subprocess.run(
        [
            "snakemake",
            "--cores",
            "1",
            "--verbose",
            "--printshellcmds",
            "--keep-going",
            "-s",
            str(WORKFLOW_PATH),
            "--configfile",
            str(snake_test_env["configfile"]),
        ],
        cwd=snake_test_env["workdir"],
        capture_output=True,
        text=True,
    )

    print(result.stdout)
    assert (
        result.returncode == 0
    ), f"Snakemake execution failed:\n{result.stderr}"

    # Check that output files were created
    output_base = output_base = (
        snake_test_env["workdir"]
        / "processed"
        / "sub-0_001"
        / "ses-0"
        / "funcimg"
        / "suite2p"
        / "plane0"
    )
    assert (output_base / "stat.npy").exists(), "Missing output: stat.npy"
    assert (output_base / "data.bin").exists(), "Missing output: data.bin"
