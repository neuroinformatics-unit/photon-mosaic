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
        "use_slurm": False,
        "suite2p_ops": {
            "fs": 6.0,
            "nplanes": 1,
            "tau": 0.7,
            "nonrigid": True,
            "diameter": 8,
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
