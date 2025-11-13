"""
Unit tests for the CLI module, specifically testing the automatic
unlock and retry functionality.
"""

import platform
import signal
import subprocess
import time


def test_automatic_unlock_on_interrupted_workflow(snake_test_env):
    """Test that photon-mosaic automatically handles locked workflows.

    This test verifies the real-world scenario where:
    1. A photon-mosaic process is interrupted, leaving the directory locked
    2. Re-running photon-mosaic automatically detects the lock
    3. The lock is automatically removed
    4. The pipeline successfully executes

    This simulates what happens when a user Ctrl+C's a running workflow
    or when a process is killed unexpectedly.
    """
    workdir = snake_test_env["workdir"]
    configfile = snake_test_env["configfile"]

    # Build the photon-mosaic command
    cmd = [
        "photon-mosaic",
        "--config",
        str(configfile),
        "--jobs",
        "1",
    ]

    print("\n=== Step 1: Starting photon-mosaic and interrupting it ===")
    print(f"Command: {' '.join(cmd)}")
    print(f"Working directory: {workdir}")

    # Start the first process
    # On Windows, we need to create a new process group to avoid
    # signals affecting the parent process
    creation_flags = 0
    if platform.system() == "Windows":
        creation_flags = subprocess.CREATE_NEW_PROCESS_GROUP

    with subprocess.Popen(
        cmd,
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        creationflags=creation_flags,
    ) as process:
        # Give it time to start and create locks
        # (Snakemake creates locks early in execution)
        time.sleep(2)

        # Interrupt the process to leave it in a locked state
        print("Interrupting process to create a locked state...")

        # On Windows, CTRL_BREAK_EVENT would affect the parent process too,
        # so we use terminate() instead to cleanly stop the subprocess
        # On Unix, we can use SIGINT to simulate Ctrl+C
        if platform.system() == "Windows":
            # terminate() sends SIGTERM on Windows, which is cleaner
            # than CTRL_BREAK_EVENT and doesn't affect parent
            process.terminate()
        else:
            process.send_signal(signal.SIGINT)

        # Wait for process to terminate
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            # Force kill if it doesn't respond to the interrupt signal
            process.kill()
            process.wait()

        print(f"Process terminated with return code: {process.returncode}")

    # Verify that the .snakemake directory exists
    # (indicating snakemake started)
    snakemake_dir = workdir / ".snakemake"
    print(f"\nChecking for .snakemake directory: {snakemake_dir}")
    print(f".snakemake exists: {snakemake_dir.exists()}")

    if snakemake_dir.exists():
        print("Contents of .snakemake directory:")
        for item in snakemake_dir.rglob("*"):
            print(f"  {item.relative_to(snakemake_dir)}")

    # Check for lock-related files
    locks_dir = snakemake_dir / "locks"
    incomplete_dir = snakemake_dir / "incomplete"
    print(f"\nLocks directory exists: {locks_dir.exists()}")
    print(f"Incomplete directory exists: {incomplete_dir.exists()}")

    print("\n=== Step 2: Re-running photon-mosaic (should auto-unlock) ===")

    # Run the command again - it should automatically unlock and succeed
    result = subprocess.run(
        cmd,
        cwd=workdir,
        capture_output=True,
        text=True,
    )

    print(f"\nSecond run return code: {result.returncode}")
    print("\n--- STDOUT ---")
    print(result.stdout)
    print("\n--- STDERR ---")
    print(result.stderr)

    # Check the log file for evidence of auto-unlock
    logs_dir = workdir / "derivatives" / "photon-mosaic" / "logs"
    print(f"\n=== Step 3: Checking logs in {logs_dir} ===")

    if logs_dir.exists():
        log_files = list(logs_dir.glob("*.log"))
        print(f"Found {len(log_files)} log file(s)")

        for log_file in sorted(log_files):
            print(f"\nReading log file: {log_file.name}")
            log_content = log_file.read_text()

            # Print relevant sections
            if "lock" in log_content.lower():
                print("Found lock-related content in log:")
                for line in log_content.split("\n"):
                    if "lock" in line.lower():
                        print(f"  {line}")

            if "Auto-unlock attempt" in log_content:
                print("Found auto-unlock marker in log")

            if "Retry after unlock" in log_content:
                print("Found retry marker in log")

        # Now verify the second run succeeded
        assert result.returncode == 0, (
            f"photon-mosaic failed after attempted auto-unlock.\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}\n"
            f"Check logs in: {logs_dir}"
        )

        # Verify that the log contains evidence of automatic unlock
        # (only if the workflow was actually locked)
        latest_log = max(log_files, key=lambda p: p.stat().st_mtime)
        log_content = latest_log.read_text()

        # Check if there was a lock error that triggered auto-unlock
        has_lock_error = "cannot be locked" in log_content.lower()
        has_auto_unlock = "Auto-unlock attempt" in log_content
        has_retry = "Retry after unlock" in log_content

        if has_lock_error:
            # If there was a lock error, verify auto-unlock was attempted
            assert (
                has_auto_unlock
            ), "Lock error detected but no auto-unlock attempt found in log"
            assert (
                has_retry
            ), "Lock error detected but no retry marker found in log"
            print(
                "\n✓ Verified: Lock was detected and "
                "automatic unlock was performed"
            )
        else:
            print(
                "\n✓ No lock detected (process may not have "
                "held lock long enough)"
            )

        print("\n✓ Test passed: Pipeline completed successfully")
    else:
        # If logs don't exist, that's also acceptable for a successful run
        assert result.returncode == 0, (
            f"photon-mosaic failed.\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
        print("\n✓ Test passed: Pipeline completed successfully")
