"""
Regression tests for file structure consistency.

This module compares the file structure logs generated during test runs
against saved expected logs to catch unexpected changes in file structure
generation behavior.
"""

import re
from pathlib import Path

import pytest


@pytest.mark.order(-1)  # Run this test last
def test_file_structure_regression():
    """
    Compare current file structure logs against expected baseline.

    This test runs last (after all other tests have generated their logs)
    and compares the file structures against saved expected logs. If the
    expected logs don't exist yet, they are copied from the current logs.
    """
    test_root = Path(__file__).parent.parent
    logs_dir = test_root / "logs"
    expected_logs_dir = test_root / "expected_logs"

    # Ensure directories exist
    logs_dir.mkdir(exist_ok=True)
    expected_logs_dir.mkdir(exist_ok=True)

    # Get all current log files
    current_logs = list(logs_dir.glob("*.log"))

    # Group logs by test name (remove timestamp suffix)
    test_logs = {}
    for log_file in current_logs:
        # Extract test name (everything before the timestamp)
        # Timestamp format is YYYYMMDD_HHMMSS_microseconds (3 parts)
        test_name = "_".join(log_file.stem.split("_")[:-3])
        if test_name not in test_logs:
            test_logs[test_name] = []
        test_logs[test_name].append(log_file)

    mismatches = []
    newly_created = []

    def normalize_log_content(log_path):
        """
        Read log and normalize for comparison.

        Removes:
        - Timestamp lines (vary between runs)
        - Root path lines (vary between systems/runs)
        - Timestamps in filenames (YYYYMMDD_HHMMSS format)
        - Full absolute paths (replace with relative paths)
        """
        with open(log_path, "r") as f:
            lines = f.readlines()

        normalized = []
        timestamp_pattern = re.compile(r"_\d{8}_\d{6}")
        # Pattern to match absolute paths (both Unix and Windows style)
        abs_path_pattern = re.compile(r"(/[^\s:]+|[A-Z]:\\[^\s:]+)")

        for line in lines:
            # Skip timestamp and root path lines entirely
            if line.startswith("Timestamp:") or line.startswith("Root:"):
                continue

            # Replace timestamp patterns in filenames with placeholder
            normalized_line = timestamp_pattern.sub("_TIMESTAMP", line)

            # Replace absolute paths with placeholder
            # This handles paths in any context (not just Root: lines)
            normalized_line = abs_path_pattern.sub("<PATH>", normalized_line)

            normalized.append(normalized_line)

        return "".join(normalized)

    for test_name, log_files in test_logs.items():
        # Use the most recent log for this test
        latest_log = sorted(log_files, reverse=True)[0]
        expected_log = expected_logs_dir / f"{test_name}.log"

        if not expected_log.exists():
            # First run - save normalized version as expected log
            normalized_content = normalize_log_content(latest_log)
            with open(expected_log, "w") as f:
                f.write(normalized_content)
            newly_created.append(test_name)
        else:
            # Compare file structure (ignore timestamps and root paths)
            current_content = normalize_log_content(latest_log)
            expected_content = normalize_log_content(expected_log)
            if current_content != expected_content:
                mismatches.append(test_name)

    if mismatches:
        error_msg = (
            f"\n❌ Found {len(mismatches)} file structure mismatch(es):\n\n"
        )
        for test_name in mismatches:
            error_msg += f"  - {test_name}\n"

        error_msg += "\nTo update expected logs after intentional changes:\n"
        error_msg += f"  rm -rf {expected_logs_dir}/*.log\n"
        error_msg += "  pytest\n"

        pytest.fail(error_msg)
