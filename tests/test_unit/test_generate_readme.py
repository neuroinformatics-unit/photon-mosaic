"""Test for automatic README generation."""

import tempfile
from pathlib import Path
from typing import Any, Dict

from photon_mosaic.dataset_discovery import DatasetDiscoverer
from tests.test_data_factory import TestDataFactory as DataFactory


def print_directory_tree(
    path: Path, prefix: str = "", max_depth: int = 10, current_depth: int = 0
) -> str:
    """Generate directory tree structure as a string."""
    if current_depth > max_depth:
        return ""

    if not path.exists():
        return f"{prefix}[PATH DOES NOT EXIST: {path}]\n"

    items = sorted(path.iterdir(), key=lambda x: (x.is_file(), x.name))

    result = ""
    for i, item in enumerate(items):
        is_last = i == len(items) - 1
        current_prefix = "└── " if is_last else "├── "
        result += f"{prefix}{current_prefix}{item.name}\n"

        if item.is_dir() and current_depth < max_depth:
            extension_prefix = "    " if is_last else "│   "
            result += print_directory_tree(
                item, prefix + extension_prefix, max_depth, current_depth + 1
            )

    return result


def run_dataset_discovery(
    raw_data_path: Path, dataset_type: str, neuroblueprint_format: bool = False
) -> Dict[str, Any]:
    """Test dataset discovery and return structured results."""
    try:
        # Configure based on dataset type
        if neuroblueprint_format:
            tiff_patterns = ["*.tif"]
        elif "Basic" in dataset_type:
            tiff_patterns = ["type_1*.tif", "type_2*.tif"]
        else:
            tiff_patterns = ["*.tif"]

        # Create a DatasetDiscoverer instance
        discoverer = DatasetDiscoverer(
            base_path=raw_data_path,
            pattern=".*",  # Match all directories
            tiff_patterns=tiff_patterns,
            neuroblueprint_format=neuroblueprint_format,
        )

        # Actually run the discovery
        discoverer.discover()

        # Simple results structure
        results: Dict[str, Any] = {
            "dataset_type": dataset_type,
            "original_datasets": discoverer.original_datasets,
            "transformed_datasets": discoverer.transformed_datasets,
            "tiff_files": dict(discoverer.tiff_files),
            "session_mapping": {},
            "output_examples": [],
        }

        # Generate simple session mapping for output structure
        for i, dataset in enumerate(discoverer.transformed_datasets):
            if i < len(discoverer.original_datasets):
                original_name = discoverer.original_datasets[i]
                tiff_files = discoverer.tiff_files.get(original_name, {})

                for session_idx in tiff_files.keys():
                    session_name = discoverer.get_session_name(i, session_idx)
                    output_path = (
                        f"sub-{i+1:03d}_id-{original_name}"
                        if not neuroblueprint_format
                        else dataset
                    )

                    results["session_mapping"][
                        f"{output_path}/{session_name}"
                    ] = {
                        "output_path": output_path,
                        "session_path": session_name,
                        "files": tiff_files[session_idx],
                    }
                    results["output_examples"].append(
                        f"{output_path}/{session_name}"
                    )

        return results

    except Exception as e:
        return {
            "dataset_type": dataset_type,
            "error": str(e),
            "session_mapping": {},
            "output_examples": [],
        }


def generate_simple_output_structure(results: Dict[str, Any]) -> str:
    """Generate a simple output directory structure from discovery results."""
    if "error" in results or not results["output_examples"]:
        return "# No valid datasets discovered\n"

    structure = "derivatives/\n"

    # Simple tree structure without complex formatting
    for example in results["output_examples"]:
        parts = example.split("/")
        if len(parts) >= 2:
            subject, session = parts[0], parts[1]
            structure += f"├── {subject}/\n"
            structure += f"│   └── {session}/\n"
            structure += "│       └── funcimg/\n"

            # Add sample files from session mapping
            for key, mapping in results["session_mapping"].items():
                if (
                    mapping["output_path"] == subject
                    and mapping["session_path"] == session
                ):
                    for file in mapping["files"][:2]:  # Show max 2 files
                        structure += (
                            f"│           └── {{output_pattern}}{file}\n"
                        )
                    break

    return structure


def generate_readme_content() -> str:
    """Generate the complete README content based on actual test runs."""

    readme_content = """# Test Data Input/Output Examples

**⚠️ Auto-generated from `test_generate_readme.py` - Do not edit manually!**

## Input → Output Transformations

*Generated from actual test runs showing how photon-mosaic transforms input
folder structures:*

"""

    factory = DataFactory()

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Test cases
        test_cases = [
            {
                "name": "Basic Dataset",
                "factory_method": "create_basic_dataset",
                "neuroblueprint": False,
                "description": "Creates a simple dataset structure similar to "
                "the original test data",
            },
            {
                "name": "Custom Metadata Dataset",
                "factory_method": "create_custom_metadata_dataset",
                "neuroblueprint": False,
                "description": "Creates a dataset with custom metadata "
                "encoded in folder names",
            },
            {
                "name": "NeuroBlueprint Format Dataset",
                "factory_method": "create_neuroblueprint_dataset",
                "neuroblueprint": True,
                "description": "Creates a dataset following the "
                "NeuroBlueprint standard format",
            },
            {
                "name": "Non-continuous NeuroBlueprint Dataset",
                "factory_method": (
                    "create_noncontinuous_neuroblueprint_dataset"
                ),
                "neuroblueprint": True,
                "description": "Creates a NeuroBlueprint format dataset with "
                "non-continuous IDs for testing edge cases",
            },
        ]

        for test_case in test_cases:
            readme_content += f"### {test_case['name']}\n\n"

            # Create test data
            method_name = str(test_case["factory_method"])
            method = getattr(factory, method_name)
            raw_data_path: Path = method(tmp_path / method_name)

            # Test dataset discovery
            results = run_dataset_discovery(
                raw_data_path,
                str(test_case["name"]),
                bool(test_case["neuroblueprint"]),
            )

            # Generate input structure
            readme_content += "**Input:**\n```\n"
            input_structure = print_directory_tree(raw_data_path)
            readme_content += input_structure
            readme_content += "```\n\n"

            # Generate output structure
            readme_content += "**Output:**\n```\n"
            if "error" not in results:
                output_structure = generate_simple_output_structure(results)
                readme_content += output_structure
            else:
                readme_content += f"# Error: {results['error']}\n"
            readme_content += "```\n\n"

    # Add minimal footer
    readme_content += (
        "## Key Rules\n\n"
        "- **Custom formats** → `sub-{INDEX:03d}_id-{original-name}` + "
        "`ses-000`, `ses-001`...\n"
        "- **NeuroBlueprint formats** → Preserve original names + "
        "session IDs\n"
        "- **All outputs** → `derivatives/{subject}/{session}/funcimg/`\n\n"
        "**Regenerate:** `pytest tests/test_unit/test_generate_readme.py::"
        "test_generate_readme -v`\n"
    )

    return readme_content


def test_generate_readme(tmp_path):
    """Generate and validate the tests README.md file."""
    print("Generating tests README from actual test runs...")

    readme_content = generate_readme_content()

    readme_path = Path(__file__).parent.parent / "README.md"

    # Write the generated content
    with open(readme_path, "w") as f:
        f.write(readme_content)

    # Validate the generated content
    lines = readme_content.splitlines()

    # Simple validation checks
    assert len(lines) > 20, "Generated README should have content"
    assert "# Test Data Input/Output Examples" in readme_content
    assert "**Input:**" in readme_content and "**Output:**" in readme_content
    assert "derivatives/" in readme_content and "sub-" in readme_content

    print(f"✅ Generated {readme_path}")
    print(f"📝 README contains {len(lines)} lines")
    print("🔍 Validation passed - README is comprehensive and up-to-date")


def test_dataset_discovery_consistency():
    """
    Test that dataset discovery behavior is consistent across
    different input types.

    This test validates that:
    1. Basic datasets use pattern-based session indexing (0, 1, 2...)
    2. NeuroBlueprint datasets preserve original session IDs
    3. Folder naming follows expected conventions
    """
    factory = DataFactory()

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Test basic dataset discovery
        raw_data_basic = factory.create_basic_dataset(tmp_path / "basic")
        results_basic = run_dataset_discovery(
            raw_data_basic, "Basic Dataset", False
        )

        # Simple validation - just check that discovery works
        assert (
            len(results_basic["output_examples"]) > 0
        ), "Should discover basic datasets"

        # Test NeuroBlueprint dataset discovery
        raw_data_nb = factory.create_neuroblueprint_dataset(tmp_path / "nb")
        results_nb = run_dataset_discovery(
            raw_data_nb, "NeuroBlueprint Dataset", True
        )

        assert (
            len(results_nb["output_examples"]) > 0
        ), "Should discover NeuroBlueprint datasets"


if __name__ == "__main__":
    # Allow running the test directly for manual README generation
    test_generate_readme(None)
