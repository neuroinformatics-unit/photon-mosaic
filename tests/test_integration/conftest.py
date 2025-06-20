import shutil
from pathlib import Path

import pytest
import yaml


@pytest.fixture
def base_config():
    """Create a base configuration that can be extended."""

    photon_mosaic_path = Path(__file__).parent.parent.parent
    with open(
        photon_mosaic_path / "photon_mosaic" / "workflow" / "config.yaml", "r"
    ) as f:
        config = yaml.safe_load(f)

    config["raw_data_base"] = "raw_data"
    config["processed_data_base"] = "derivatives"
    #  match different file naming patterns to different sessions
    config["dataset_discovery"]["tiff_patterns"] = [
        "type_1*.tif",
        "type_2*.tif",
    ]

    return config


@pytest.fixture
def map_of_tiffs():
    """Create a map of tiffs in the test data using rglob"""
    photon_mosaic_path = Path(__file__).parent.parent.parent / "tests" / "data"
    map_of_tiffs = {}
    for dataset in photon_mosaic_path.glob("*"):
        if dataset.is_dir():
            # Get just the filenames, not the full paths
            tiff_files = [f.name for f in dataset.rglob("*.tif")]
            map_of_tiffs[dataset.name] = tiff_files
    return map_of_tiffs


@pytest.fixture
def test_config_with_contrast(base_config):
    """
    Create a test configuration with contrast enhancement preprocessing step.
    """
    config = base_config.copy()
    config["preprocessing"] = {
        "steps": [
            {
                "name": "contrast",
                "kwargs": {
                    "percentile_low": 1.0,
                    "percentile_high": 99.0,
                },
            }
        ],
        "output_pattern": "enhanced_",
    }
    return config


@pytest.fixture
def snake_test_env(tmp_path, base_config):
    """
    Fixture that sets up the test environment with data and configuration.
    """
    print("\n=== Setting up test environment ===")
    print(f"Temporary directory: {tmp_path}")

    raw_data = tmp_path / "raw_data"
    print(f"Raw data directory: {raw_data}")
    print(f"Copying test data from: {Path('tests/data').absolute()}")
    shutil.copytree("tests/data", raw_data)
    print(f"Raw data contents after copy: {list(raw_data.glob('**/*'))}")

    processed_data = tmp_path / "derivatives"
    processed_data.mkdir()
    print(f"Processed data directory: {processed_data}")

    # Update paths in config
    config = base_config.copy()
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
