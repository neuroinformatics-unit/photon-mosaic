"""
Shared fixtures for all tests in the photon-mosaic test suite.

This module provides common fixtures used across both unit and integration
tests, following the DRY principle to avoid duplication.
"""

from pathlib import Path

import pytest
import yaml

from .test_data_factory import TestDataFactory


@pytest.fixture
def test_data_root():
    """Return the path to test data directory."""
    return Path(__file__).parent


@pytest.fixture
def data_factory():
    """Return a TestDataFactory instance for creating test data dynamically."""
    return TestDataFactory()


@pytest.fixture
def base_config():
    """Create a base configuration that can be extended."""
    photon_mosaic_path = Path(__file__).parent.parent
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
def metadata_base_config():
    """Create a base configuration for metadata testing."""
    photon_mosaic_path = Path(__file__).parent.parent
    with open(
        photon_mosaic_path / "photon_mosaic" / "workflow" / "config.yaml", "r"
    ) as f:
        config = yaml.safe_load(f)

    # Set default values for metadata testing
    config["dataset_discovery"]["tiff_patterns"] = ["*.tif"]
    return config


@pytest.fixture
def map_of_tiffs():
    """
    Create a map of tiffs in test data using rglob -
    for backward compatibility
    """

    photon_mosaic_path = Path(__file__).parent / "data"
    map_of_tiffs = {}
    for dataset in photon_mosaic_path.glob("*"):
        if dataset.is_dir():
            # Get just the filenames, not the full paths
            tiff_files = [f.name for f in dataset.rglob("*.tif")]
            map_of_tiffs[dataset.name] = tiff_files
    return map_of_tiffs


def create_map_of_tiffs(raw_data_path: Path) -> dict:
    """
    Create a map of tiffs for a given raw data directory.

    Args:
        raw_data_path: Path to raw data directory

    Returns:
        Dictionary mapping dataset names to list of TIFF filenames
    """
    map_of_tiffs = {}
    for dataset in raw_data_path.glob("*"):
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
def snake_test_env(tmp_path, base_config, data_factory):
    """
    Fixture that sets up the test environment with data and configuration.
    """
    print("\n=== Setting up test environment ===")
    print(f"Temporary directory: {tmp_path}")

    # Use factory to create basic dataset structure dynamically
    raw_data = data_factory.create_basic_dataset(tmp_path)
    print(f"Raw data directory: {raw_data}")
    print(f"Raw data contents after creation: {list(raw_data.glob('**/*'))}")

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


@pytest.fixture
def custom_metadata_env(tmp_path, metadata_base_config, data_factory):
    """Set up test environment for custom metadata format."""
    # Use factory to create custom metadata dataset dynamically
    raw_data = data_factory.create_custom_metadata_dataset(tmp_path)

    processed_data = tmp_path / "derivatives"
    processed_data.mkdir()

    # Update config for custom format
    config = metadata_base_config.copy()
    config["raw_data_base"] = str(raw_data.resolve())
    config["processed_data_base"] = str(processed_data.resolve())
    config["dataset_discovery"]["neuroblueprint_format"] = False
    config["dataset_discovery"]["pattern"] = ".*"

    # Create config file
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f)

    return {
        "workdir": tmp_path,
        "configfile": config_path,
        "raw_data": raw_data,
        "processed_data": processed_data,
    }


@pytest.fixture
def neuroblueprint_env(tmp_path, metadata_base_config, data_factory):
    """Set up test environment for NeuroBlueprint metadata format."""
    # Use factory to create NeuroBlueprint dataset dynamically
    raw_data = data_factory.create_neuroblueprint_dataset(tmp_path)

    processed_data = tmp_path / "derivatives"
    processed_data.mkdir()

    # Update config for NeuroBlueprint format
    config = metadata_base_config.copy()
    config["raw_data_base"] = str(raw_data.resolve())
    config["processed_data_base"] = str(processed_data.resolve())
    config["dataset_discovery"]["neuroblueprint_format"] = True

    # Create config file
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f)

    return {
        "workdir": tmp_path,
        "configfile": config_path,
        "raw_data": raw_data,
        "processed_data": processed_data,
    }


@pytest.fixture
def neuroblueprint_noncontinuous_env(
    tmp_path, metadata_base_config, data_factory
):
    """
    Set up test environment for NeuroBlueprint format with non-continuous IDs.
    """
    # Use factory to create non-continuous ID dataset dynamically
    raw_data = data_factory.create_noncontinuous_neuroblueprint_dataset(
        tmp_path
    )

    processed_data = tmp_path / "derivatives"
    processed_data.mkdir()

    # Update config for NeuroBlueprint format
    config = metadata_base_config.copy()
    config["raw_data_base"] = str(raw_data.resolve())
    config["processed_data_base"] = str(processed_data.resolve())
    config["dataset_discovery"]["neuroblueprint_format"] = True

    # Create config file
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f)

    return {
        "workdir": tmp_path,
        "configfile": config_path,
        "raw_data": raw_data,
        "processed_data": processed_data,
    }
