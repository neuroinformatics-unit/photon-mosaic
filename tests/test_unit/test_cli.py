"""
Unit tests for the CLI module, specifically testing the
execute_pipeline_with_retry function.
"""

from unittest.mock import MagicMock, patch

import pytest

from photon_mosaic.cli import execute_pipeline_with_retry


@pytest.fixture
def temp_log_path(tmp_path):
    """Create a temporary log file path."""
    return tmp_path / "test.log"


@pytest.fixture
def mock_snakemake_workspace(tmp_path):
    """Create a temporary directory structure mimicking Snakemake's
    workspace with lock directory."""
    workspace = tmp_path / "snakemake_workspace"
    workspace.mkdir()
    snakemake_dir = workspace / ".snakemake"
    snakemake_dir.mkdir()
    locks_dir = snakemake_dir / "locks"
    locks_dir.mkdir()
    return workspace, locks_dir


def test_execute_pipeline_with_lock_and_retry(
    temp_log_path, mock_snakemake_workspace
):
    """Test that photon-mosaic automatically handles locked workflows.

    This test verifies the full lock detection and retry cycle:
    1. Initial command fails with lock error
    2. Function detects lock in output
    3. Unlock command is executed
    4. Original command is retried and succeeds
    """
    workspace, locks_dir = mock_snakemake_workspace
    lock_file = locks_dir / "test.lock"

    call_count = {"count": 0}

    def mock_subprocess_run(cmd, stdout=None, stderr=None, **kwargs):
        """Mock subprocess that simulates realistic Snakemake behavior."""
        call_count["count"] += 1

        # First call: workflow fails due to existing lock
        if call_count["count"] == 1:
            # Simulate Snakemake creating a lock file
            lock_file.write_text("locked by process 12345")

            # Write realistic Snakemake lock error to log
            if stdout:
                stdout.write(
                    "Error: Directory cannot be locked. Please make "
                    "sure that no other Snakemake process is running "
                    "on this directory.\n"
                )
            result = MagicMock(returncode=1)
            return result

        # Second call: unlock command
        elif call_count["count"] == 2:
            # Verify this is the unlock command
            assert "--unlock" in cmd, "Expected unlock command"

            # Simulate successful unlock by removing lock file
            if lock_file.exists():
                lock_file.unlink()

            if stdout:
                stdout.write("Unlocking working directory.\n")
            result = MagicMock(returncode=0)
            return result

        # Third call: retry succeeds (lock is gone)
        else:
            # Verify unlock actually removed the lock
            assert not lock_file.exists(), "Lock file should be removed"

            if stdout:
                stdout.write("Workflow completed successfully\n")
            result = MagicMock(returncode=0)
            return result

    with patch("photon_mosaic.cli.subprocess.run") as mock_run:
        mock_run.side_effect = mock_subprocess_run

        with patch("photon_mosaic.cli.logging.getLogger") as mock_get_logger:
            mock_logger = MagicMock()
            mock_get_logger.return_value = mock_logger

            cmd = ["snakemake", "--snakefile", "test.smk"]
            result = execute_pipeline_with_retry(cmd, temp_log_path)

            # Verify successful completion after automatic unlock
            assert result == 0
            assert call_count["count"] == 3

            # Verify the correct sequence of commands was executed
            calls = mock_run.call_args_list
            assert len(calls) == 3

            # First call: original command
            assert calls[0][0][0] == cmd

            # Second call: unlock command
            assert calls[1][0][0] == cmd + ["--unlock"]

            # Third call: retry original command
            assert calls[2][0][0] == cmd

            # Verify appropriate log messages
            mock_logger.warning.assert_called_once()
            warning_msg = mock_logger.warning.call_args[0][0]
            assert "locked" in warning_msg.lower()
            assert "automatic" in warning_msg.lower()

            mock_logger.info.assert_any_call(
                "Automatic unlock successful. Retrying pipeline execution..."
            )

            # Verify log file structure
            with open(temp_log_path, "r") as f:
                log_content = f.read()

            assert "cannot be locked" in log_content.lower()
            assert "--- Auto-unlock attempt ---" in log_content
            assert "Unlocking working directory" in log_content
            assert "--- Retry after unlock ---" in log_content
            assert "completed successfully" in log_content.lower()
