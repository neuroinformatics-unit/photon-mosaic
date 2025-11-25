# File Structure Regression Tests

This directory contains regression tests that verify the consistency of file structures generated during test runs.

## How It Works

1. **Log Generation**: As tests run, the `log_test_fs` fixture (in `conftest.py`) automatically captures the file structure after each test and saves it to `tests/logs/`.

2. **Expected Logs**: On the first run, the regression test (`test_file_structure_regression.py`) copies the logs to `tests/expected_logs/` as the baseline.

3. **Comparison**: On subsequent runs, the regression test compares new logs against the expected logs and fails if there are differences.

## Running the Tests

```bash
# Run all tests (including regression test at the end)
pytest

# Run only the regression test
pytest tests/test_regression/test_file_structure_regression.py
```

## Updating Expected Logs

If you've intentionally changed the file structure generation behavior, update the expected logs:

```bash
# Remove old expected logs
rm -rf tests/expected_logs/*.log

# Run tests to generate new expected logs
pytest
```
