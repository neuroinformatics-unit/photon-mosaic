"""
Tests for metadata functionality in dataset discovery.

This module tests the metadata extraction capabilities of the DatasetDiscoverer
for both custom metadata and NeuroBlueprint formats.
"""

import subprocess

from photon_mosaic.dataset_discovery import DatasetDiscoverer


def run_photon_mosaic_dry_run(workdir, configfile):
    """Helper function to run photon-mosaic CLI with dry-run."""
    cmd = [
        "photon-mosaic",
        "--config",
        str(configfile),
        "--dry-run",
        "--log-level",
        "DEBUG",
    ]

    result = subprocess.run(
        cmd, cwd=workdir, capture_output=True, text=True, timeout=60
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

        # Check that we found the expected dataset
        original_names = discoverer.original_datasets
        assert (
            "mouse001_genotype-WT_age-P60_treatment-saline" in original_names
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

        # Check that we found expected datasets
        original_names = discoverer.original_datasets
        expected_names = [
            "sub-001_strain-C57BL6_sex-M",
            "sub-002_strain-C57BL6_sex-F",
        ]
        assert any(name in original_names for name in expected_names)

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
        result = run_photon_mosaic_dry_run(
            custom_metadata_env["workdir"], custom_metadata_env["configfile"]
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
        result = run_photon_mosaic_dry_run(
            neuroblueprint_env["workdir"], neuroblueprint_env["configfile"]
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
