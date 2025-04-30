(user_guide/dataset_discovery)=
# Dataset Discovery

The `photon-mosaic` pipeline includes a flexible dataset discovery system that allows you to automatically find and process datasets based on configurable patterns and transformations.

## Overview

The dataset discovery system uses regular expressions (regex) to:
1. Find datasets matching a specific pattern
2. Exclude datasets matching certain patterns
3. Transform dataset names using regex substitutions

This makes it easy to adapt the pipeline to different dataset naming conventions and organizational structures.

## Configuration

Configure dataset discovery in your `config.yaml` file under the `dataset_discovery` section:

```yaml
dataset_discovery:
  # Regex pattern to match dataset names
  pattern: "^2.*"  # Example: match directories starting with '2'

  # Regex patterns for datasets to exclude
  exclude_patterns:
    - ".*_test$"  # Exclude test datasets
    - ".*_backup$"  # Exclude backup datasets

  # Transform names using regex substitutions
  substitutions:
    # Remove underscores
    - pattern: "_"
      repl: ""

    # Take the first part of the path (e.g., "dataset/subfolder" -> "dataset")
    - pattern: "([^/]+)/.*"
      repl: "\1"
```

## Pattern Matching

The `pattern` field uses a regular expression to match dataset names. Some examples:

- `^2.*` - Match directories starting with '2' (e.g., date-based directories like "2023_01")
- `.*_experiment_.*` - Match directories containing "_experiment_"
- `exp_\d+` - Match directories like "exp_001", "exp_002", etc.

## Excluding Datasets

The `exclude_patterns` field allows you to exclude datasets matching certain patterns:

- `.*_test$` - Exclude directories ending with "_test"
- `.*_backup$` - Exclude directories ending with "_backup"
- `temp_.*` - Exclude directories starting with "temp_"

## Transforming Dataset Names

The `substitutions` field allows you to transform dataset names using regex substitutions. Each substitution has:

- `pattern`: The regex pattern to match
- `repl`: The replacement string (can use capture groups with `\1`, `\2`, etc.)

Common transformations:

```yaml
# Remove underscores
- pattern: "_"
  repl: ""

# Convert date format from "2023_01" to "202301"
- pattern: "(\d{4})_(\d{2})"
  repl: "\1\2"

# Take only the first part of a path
- pattern: "([^/]+)/.*"
  repl: "\1"
```

## Advanced Examples

### Date-based Datasets

For datasets organized by date:

```yaml
dataset_discovery:
  pattern: "^2\d{3}_\d{2}_\d{2}"  # Match YYYY_MM_DD format
  exclude_patterns:
    - ".*_test$"
  substitutions:
    - pattern: "_"
      repl: ""
```

### Experiment-based Datasets

For datasets organized by experiment:

```yaml
dataset_discovery:
  pattern: "exp_\d+"
  exclude_patterns:
    - ".*_backup$"
  substitutions:
    - pattern: "exp_(\d+)"
      repl: "experiment\1"
```

### Complex Transformations

For more complex naming conventions:

```yaml
dataset_discovery:
  pattern: ".*"
  exclude_patterns:
    - "temp_.*"
    - ".*_old$"
  substitutions:
    # Remove version suffixes
    - pattern: "_v\d+$"
      repl: ""
    # Convert date format
    - pattern: "(\d{4})_(\d{2})_(\d{2})"
      repl: "\1\2\3"
    # Remove specific prefixes
    - pattern: "^prefix_"
      repl: ""
```

## Integration with Snakemake

The dataset discovery system is integrated with Snakemake and automatically used in the workflow. The discovered datasets are available as wildcards in your rules:

```python
rule all:
    input:
        expand(
            "results/{dataset}/processed.tif",
            dataset=datasets  # 'datasets' is automatically populated
        )
```

You don't need to modify your Snakefile to use the dataset discovery system - it's automatically configured based on your `config.yaml` file.
