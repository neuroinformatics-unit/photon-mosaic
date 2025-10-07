"""
Tests for metadata functionality in dataset discovery.

This module tests the extraction capabilities of the DatasetDiscoverer
for both custom metadata and NeuroBlueprint formats.
"""

import re
import shutil
import subprocess
from pathlib import Path

import yaml

from photon_mosaic.dataset_discovery import DatasetDiscoverer


def run_photon_mosaic(workdir, configfile, timeout=None):
    """Helper function to run photon-mosaic CLI with dry-run.

    timeout: seconds to wait for the subprocess to complete. If None,
    wait indefinitely (no timeout).
    """
    cmd = [
        "photon-mosaic",
        "--config",
        str(configfile),
        "--log-level",
        "DEBUG",
    ]

    result = subprocess.run(
        cmd, cwd=workdir, capture_output=True, text=True, timeout=timeout
    )

    return result


class TestMetadataFunctionality:
    """Test class for metadata extraction functionality."""

    def test_custom_metadata_extraction(self, custom_metadata_env):
        """Test metadata extraction from custom format data."""
        # Create discoverer for custom metadata format
        discoverer = DatasetDiscoverer(
            base_path=custom_metadata_env["raw_data"],
            pattern=".*",
            tiff_patterns=["*.tif"],
            neuroblueprint_format=False,
        )

        # Discover datasets
        discoverer.discover()

        # Basic checks
        assert len(discoverer.datasets) > 0, "Should find at least one dataset"

        # Check that we found datasets with the expected metadata pattern
        original_names = discoverer.original_datasets
        # Factory creates "mouse001_genotype-WT_age-P60_treatment-saline"
        # by default
        assert any(
            "mouse" in name and "genotype-" in name for name in original_names
        ), (
            f"Expected subjects with mouse "
            f"and genotype metadata in {original_names}"
        )

        # Verify transformed names follow expected pattern
        transformed_names = discoverer.transformed_datasets
        assert len(transformed_names) == len(original_names)
        assert all(name.startswith("sub-") for name in transformed_names)

        # Check that dataset has empty subject metadata for custom format
        dataset = discoverer.datasets[0]
        assert (
            dataset.subject_metadata == ""
        ), "Custom format should have empty subject metadata"

    def test_neuroblueprint_metadata_extraction(self, neuroblueprint_env):
        """Test metadata extraction from NeuroBlueprint format data."""
        # Create discoverer for NeuroBlueprint format
        discoverer = DatasetDiscoverer(
            base_path=neuroblueprint_env["raw_data"],
            tiff_patterns=["*.tif"],
            neuroblueprint_format=True,
        )

        # Discover datasets
        discoverer.discover()

        # Basic checks
        assert (
            len(discoverer.datasets) > 0
        ), "Should find at least one NeuroBlueprint dataset"

        # Check that we found expected datasets with NeuroBlueprint format
        original_names = discoverer.original_datasets
        # The factory creates subjects with strain and sex metadata by default
        assert any(
            "sub-" in name and "strain-" in name and "sex-" in name
            for name in original_names
        ), (
            f"Expected NeuroBlueprint subjects with strain and "
            f"sex metadata in {original_names}"
        )
        # Should find at least 2 subjects (default factory behavior)
        assert (
            len(original_names) >= 2
        ), f"Expected at least 2 subjects, got {len(original_names)}"

        # Verify metadata extraction
        for dataset in discoverer.datasets:
            # Should have subject metadata extracted from folder name
            assert (
                dataset.subject_metadata != ""
            ), "NeuroBlueprint format should extract subject metadata"

            # Should contain strain and sex metadata
            subject_meta = dataset.subject_metadata
            assert (
                "strain-" in subject_meta
            ), f"Subject metadata should contain strain: {subject_meta}"
            assert (
                "sex-" in subject_meta
            ), f"Subject metadata should contain sex: {subject_meta}"

            # Check session metadata
            assert (
                len(dataset.session_metadata) > 0
            ), "Should have session metadata"

            # At least one session should have metadata
            session_metas = list(dataset.session_metadata.values())
            assert any(
                meta != "" for meta in session_metas
            ), "At least one session should have metadata"

    def test_metadata_inference(self):
        """Test that metadata keys are correctly inferred from folder names."""
        # Test the static method for inferring metadata
        folder_names = [
            "sub-001_strain-C57BL6_sex-M",
            "ses-001_date-20250225_protocol-training",
        ]

        inferred = DatasetDiscoverer._infer_metadata_keys_from_folder_names(
            folder_names
        )

        # Should find strain, sex, date, and protocol keys
        expected_keys = {"strain", "sex", "date", "protocol"}
        inferred_keys = set(inferred.keys())

        assert expected_keys.issubset(
            inferred_keys
        ), f"Expected {expected_keys}, got {inferred_keys}"

        # Check that patterns are reasonable
        assert "strain-([^_]+)" in inferred.values()
        assert "sex-([^_]+)" in inferred.values()

    def test_neuroblueprint_format_validation(self):
        """Test NeuroBlueprint format validation."""
        # Valid NeuroBlueprint subject names
        valid_subjects = [
            "sub-001",
            "sub-001_strain-C57BL6",
            "sub-001_strain-C57BL6_sex-M",
            "sub-mouse123_genotype-WT_age-P60",
        ]

        for name in valid_subjects:
            assert DatasetDiscoverer._is_neuroblueprint_format(
                name, "sub"
            ), f"{name} should be valid"

        # Valid NeuroBlueprint session names
        valid_sessions = [
            "ses-001",
            "ses-001_date-20250225",
            "ses-001_date-20250225_protocol-training",
            "ses-baseline_condition-control_paradigm-open-field",
        ]

        for name in valid_sessions:
            assert DatasetDiscoverer._is_neuroblueprint_format(
                name, "ses"
            ), f"{name} should be valid"

        # Invalid formats
        invalid_names = [
            "mouse001",  # No prefix
            "sub_001",  # Wrong separator
            "sub-",  # No identifier
            "sub-001_invalid",  # Missing value after key
            "sub-001_-value",  # Missing key before value
        ]

        for name in invalid_names:
            assert not DatasetDiscoverer._is_neuroblueprint_format(
                name, "sub"
            ), f"{name} should be invalid"

    def test_photon_mosaic_cli_custom_metadata(self, custom_metadata_env):
        """Test photon-mosaic CLI with custom metadata format."""
        # Run photon-mosaic with dry-run to test metadata processing
        result = run_photon_mosaic(
            custom_metadata_env["workdir"],
            custom_metadata_env["configfile"],
            timeout=None,
        )

        # Check that command ran successfully - this validates that the
        # metadata functionality works without crashing the pipeline
        assert result.returncode == 0, (
            f"Command failed with return code {result.returncode}. "
            f"Stderr: {result.stderr}"
        )

        # Check for successful pipeline execution
        output = result.stdout + result.stderr
        assert (
            "snakemake pipeline completed successfully" in output.lower()
        ), "Pipeline should complete successfully"

    def test_photon_mosaic_cli_neuroblueprint_metadata(
        self, neuroblueprint_env
    ):
        """Test photon-mosaic CLI with NeuroBlueprint metadata format."""
        # Run photon-mosaic with dry-run to test metadata processing
        result = run_photon_mosaic(
            neuroblueprint_env["workdir"],
            neuroblueprint_env["configfile"],
            timeout=None,
        )

        # Check that command ran successfully - this validates that the
        # NeuroBlueprint metadata functionality works without crashing the
        # pipeline
        assert result.returncode == 0, (
            f"Command failed with return code {result.returncode}. "
            f"Stderr: {result.stderr}"
        )
        # Check for successful pipeline execution
        output = result.stdout + result.stderr
        assert (
            "snakemake pipeline completed successfully" in output.lower()
        ), "Pipeline should complete successfully with NeuroBlueprint format"

    def test_noncontinuous_ids_preservation(
        self, neuroblueprint_noncontinuous_env
    ):
        """Test that non-continuous subject and session IDs are preserved."""
        # local variables use module-level imports: re, Path

        # First, discover what's actually in the test data folders
        raw_data_path = Path(neuroblueprint_noncontinuous_env["raw_data"])
        expected_data = {}

        for subject_dir in raw_data_path.iterdir():
            if subject_dir.is_dir() and subject_dir.name.startswith("sub-"):
                # Extract subject ID from folder name
                subject_match = re.match(r"sub-(\d+)", subject_dir.name)
                if subject_match:
                    # Find session folders and extract their IDs
                    session_ids = []
                    for session_dir in subject_dir.iterdir():
                        if (
                            session_dir.is_dir()
                            and session_dir.name.startswith("ses-")
                        ):
                            session_match = re.match(
                                r"ses-(\d+)", session_dir.name
                            )
                            if session_match:
                                session_ids.append(int(session_match.group(1)))

                    expected_data[subject_dir.name] = sorted(session_ids)

        # Now test discovery preserves these IDs
        # Run photon-mosaic CLI dry-run first to ensure the pipeline can be
        # invoked against this dataset (mirrors other CLI tests).
        result = run_photon_mosaic(
            neuroblueprint_noncontinuous_env["workdir"],
            neuroblueprint_noncontinuous_env["configfile"],
            timeout=None,
        )

        assert result.returncode == 0, (
            "photon-mosaic CLI dry-run failed: return code "
            f"{result.returncode}\n"
            f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

        # Continue with discovery checks
        discoverer = DatasetDiscoverer(
            base_path=neuroblueprint_noncontinuous_env["raw_data"],
            neuroblueprint_format=True,
            tiff_patterns=["*.tif"],
        )
        discoverer.discover()

        # Verify we found all expected subjects
        assert len(discoverer.datasets) == len(expected_data)

        # Verify each subject and session IDs are preserved
        for i, dataset in enumerate(discoverer.datasets):
            subject_name = dataset.original_name
            found_session_ids = sorted(dataset.tiff_files.keys())

            # Verify this subject was expected from folder scan
            assert (
                subject_name in expected_data
            ), f"Unexpected subject: {subject_name}"

            # Verify session IDs match what we found in the folders
            expected_session_ids = expected_data[subject_name]
            assert found_session_ids == expected_session_ids, (
                f"Subject {subject_name}: "
                f"expected sessions {expected_session_ids}, "
                f"got {found_session_ids}"
            )

            # Verify session names preserve original IDs (not sequential 0,1,2)
            for session_id in found_session_ids:
                session_name = discoverer.get_session_name(i, session_id)
                expected_prefix = f"ses-{session_id:03d}"
                assert session_name.startswith(expected_prefix), (
                    f"Expected session name to start with {expected_prefix}, "
                    f"got {session_name}"
                )

                # Verify it's NOT using sequential numbering
                # (ses-000, ses-001, etc)
                assert not session_name.startswith("ses-000"), (
                    f"Session name incorrectly uses "
                    f"sequential numbering: {session_name}"
                )

    def test_alphanumeric_subject_and_session_ids(
        self, tmp_path, metadata_base_config
    ):
        """Test that alphanumeric folder names are preserved in transformed
        names.

        Create a top-level subject folder with an evocative alphanumeric name
        (e.g. 'hippocampusA42') and a session folder named
        'novelEnv07' containing a TIFF in the raw_data. Run the CLI and verify
        discovery produces 'sub-001_id-hippocampusA42' and
        'ses-001_id-novelEnv07' (the transformed session name).
        """
        # use module-level imports: Path, shutil, yaml

        # Create raw data structure
        raw_data = tmp_path / "raw_data"
        # Use more evocative names for experimental mouse research
        subject_folder_name = "hippocampusA42"
        # In raw_data the session folder is just the alphanumeric name
        session_folder_name = "novelEnv07"
        session_folder = raw_data / subject_folder_name / session_folder_name
        session_folder.mkdir(parents=True, exist_ok=True)

        # Copy master tiff into the session
        master_tif = Path(__file__).parent.parent / "data" / "master.tif"
        shutil.copy2(master_tif, session_folder / "recording.tif")

        # Create derivatives folder
        processed_data = tmp_path / "derivatives"
        processed_data.mkdir()

        # Prepare config
        config = metadata_base_config.copy()
        config["raw_data_base"] = str(raw_data.resolve())
        config["processed_data_base"] = str(processed_data.resolve())
        config["dataset_discovery"]["tiff_patterns"] = ["*.tif"]
        # Keep discovery in custom mode so subject gets id-<name>
        config["dataset_discovery"]["neuroblueprint_format"] = False

        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            yaml.safe_dump(config, f)

        # Run photon-mosaic CLI (dry-run)
        result = run_photon_mosaic(tmp_path, config_path, timeout=None)
        assert result.returncode == 0, (
            f"CLI failed with return code {result.returncode}. "
            f"Stderr: {result.stderr}"
        )

        # Discover datasets and check transformed names
        discoverer = DatasetDiscoverer(
            base_path=raw_data,
            tiff_patterns=["*.tif"],
            neuroblueprint_format=False,
        )
        discoverer.discover()

        # Subject folder should be discovered as original name
        assert subject_folder_name in discoverer.original_datasets

        # Transformed subject name should include id-<subject_folder_name>
        expected_subject_transformed = f"sub-001_id-{subject_folder_name}"
        assert any(
            name == expected_subject_transformed
            for name in discoverer.transformed_datasets
        ), (
            f"Expected transformed subject name "
            f"'{expected_subject_transformed}' "
            f"got '{discoverer.transformed_datasets}'"
        )

        # Check session naming preserves the id metadata
        dataset = discoverer.datasets[0]
        found_session_ids = sorted(dataset.tiff_files.keys())
        assert found_session_ids == [1]

        session_name = discoverer.get_session_name(0, 1)
        # The discovery should transform the raw folder name into a
        # NeuroBlueprint-like session name with id metadata
        expected_session_name = "ses-001_id-novelEnv07"
        assert session_name == expected_session_name, (
            f"Expected session name '{expected_session_name}', "
            f"got '{session_name}'"
        )
