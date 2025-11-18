"""
Unit tests for the CLI module.
"""

import argparse

from photon_mosaic.cli import build_snakemake_command


def test_unlock_argument_is_appended(snake_test_env):
    """Test that --unlock is correctly appended when lock files exist.

    This test verifies that when lock files are present in the .snakemake/locks
    directory, the build_snakemake_command function correctly appends --unlock
    to the command.
    """

    workdir = snake_test_env["workdir"]
    configfile = snake_test_env["configfile"]

    # Create a fake .snakemake/locks directory to simulate a locked workflow
    snakemake_dir = workdir / ".snakemake"
    locks_dir = snakemake_dir / "locks"
    locks_dir.mkdir(parents=True, exist_ok=True)

    # Create a fake lock file to simulate a locked workflow
    lock_file = locks_dir / "0.preprocessing.lock"
    lock_file.write_text("locked")

    print(f"\n=== Created fake lock at: {lock_file} ===")

    # Create args namespace
    args = argparse.Namespace(
        config=str(configfile),
        jobs="1",
        dry_run=False,
        forcerun=None,
        rerun_incomplete=False,
        latency_wait=10,
        verbose=False,
    )

    # Build the command with lock present
    cmd_with_lock = build_snakemake_command(args, configfile, workdir)

    # Verify --unlock is in the command
    assert (
        "--unlock" in cmd_with_lock
    ), "Expected --unlock flag when lock files are present"

    # Now remove the lock file and verify --unlock is NOT added
    lock_file.unlink()
    cmd_without_lock = build_snakemake_command(args, configfile, workdir)

    # Verify --unlock is NOT in the command
    assert (
        "--unlock" not in cmd_without_lock
    ), "Did not expect --unlock flag when no lock files are present"
