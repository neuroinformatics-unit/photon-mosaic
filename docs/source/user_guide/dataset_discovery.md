(user_guide/dataset_discovery)=
# Dataset Discovery

The `photon-mosaic` pipeline includes a flexible dataset discovery system that uses regex patterns to find and process datasets.

## Configuration

Configure dataset discovery in your `config.yaml` file:

```yaml
dataset_discovery:
  # Match directories starting with '2'
  pattern: "^2.*"

  # Exclude specific patterns
  exclude_patterns:
    - ".*_test$"

  # Transform names using regex substitutions
  substitutions:
    # Remove underscores
    - pattern: "_"
      repl: ""
    # Take the first part of the path
    - pattern: "([^/]+)/.*"
      repl: "\\1"  # Note the double backslash for YAML
```

## Examples

### Date-based Datasets
```yaml
dataset_discovery:
  pattern: "^2\d{3}_\d{2}_\d{2}"  # Match YYYY_MM_DD format
  substitutions:
    - pattern: "_"
      repl: ""
```

### Experiment-based Datasets
```yaml
dataset_discovery:
  pattern: "exp_\d+"
  substitutions:
    - pattern: "exp_(\d+)"
      repl: "experiment\\1"
```

The discovered datasets are automatically used in the Snakemake workflow.
