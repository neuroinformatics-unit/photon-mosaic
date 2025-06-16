import shutil
from pathlib import Path

import pytest
import yaml


def get_base_config():
    """Create a base configuration that can be extended."""
    return {
        "raw_data_base": "raw",
        "processed_data_base": "processed",
        "dataset_discovery": {
            "pattern": "^.*$",  # Match all directories for testing
            "exclude_patterns": [],  # Don't exclude anything
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
def test_config():
    """Create a test configuration."""
    config = get_base_config()
    config["preprocessing"] = {
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
    }
    return config


@pytest.fixture
def test_config_with_contrast():
    """
    Create a test configuration with contrast enhancement preprocessing step.
    """
    config = get_base_config()
    config["preprocessing"] = {
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
    }
    return config


@pytest.fixture
def snake_test_env(tmp_path, test_config):
    """
    Fixture that sets up the test environment with data and configuration.
    """
    print("\n=== Setting up test environment ===")
    print(f"Temporary directory: {tmp_path}")

    raw_data = tmp_path / "raw"
    print(f"Raw data directory: {raw_data}")
    print(f"Copying test data from: {Path('tests/data').absolute()}")
    shutil.copytree("tests/data", raw_data)
    print(f"Raw data contents after copy: {list(raw_data.glob('**/*'))}")

    processed_data = tmp_path / "processed"
    processed_data.mkdir()
    print(f"Processed data directory: {processed_data}")

    # Update paths in config
    config = test_config.copy()
    config["raw_data_base"] = str(raw_data.resolve())
    config["processed_data_base"] = str(processed_data.resolve())

    print("\n=== Configuration ===")
    print(f"Raw data base: {config['raw_data_base']}")
    print(f"Processed data base: {config['processed_data_base']}")

    # Create config file
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, default_style='"', allow_unicode=True)
    print(f"Config file created at: {config_path}")
    print("=== End of test environment setup ===\n")

    return {
        "workdir": tmp_path,
        "configfile": config_path,
    }
