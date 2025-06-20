(user_guide/dataset_discovery)=
# Dataset Discovery

The `photon-mosaic` pipeline includes a flexible dataset discovery system that uses regex patterns to find and process datasets. This functionality allows users to process datasets with different naming conventions and structures other than the NeuroBlueprint specification.

## Configuration

Configure dataset discovery in your `config.yaml` file:

```yaml
dataset_discovery:
  # Match directories starting with '2'
  pattern: "^2.*"

  # Find tiffs. Different patterns will correspond to different sessions
  tiff_patterns: ["*.tif"]

  # Exclude specific patterns
  exclude_patterns:
    - ".*_test$"

  # Transform names using regex substitutions
  substitutions:
    - pattern: "_"
      repl: "" # The string "a_b" will be replaced with "ab"
```

## Key Parameters

- **pattern**: Regular expression to identify dataset directories
- **tiff_patterns**: List of glob patterns for TIFF files. Each pattern corresponds to a session (numbered starting from 0). Tiffs from differest sessions will be stored in different folders. Tiffs in a same session will be analysed together by Suite2p.
- **exclude_patterns**: List of regex patterns for datasets to skip during processing
- **substitutions**: List of regex substitution rules to transform dataset names

## Examples

### Single Session
```yaml
dataset_discovery:
  pattern: "^2.*"
  tiff_patterns: ["*.tif"]  # Single session with all .tif files
```

### Multiple Sessions
```yaml
dataset_discovery:
  pattern: "^2.*"
  tiff_patterns: ["session1_*.tif", "session2_*.tif"]  # Session 0, Session 1
```

### Date-based Datasets
```yaml
dataset_discovery:
  pattern: "^2\\d{3}_\\d{2}_\\d{2}"  # Match YYYY_MM_DD format
  tiff_patterns: ["*.tif"]
  substitutions:
    - pattern: "_"
      repl: ""
```

### Experiment-based Datasets
```yaml
dataset_discovery:
  pattern: "exp_\\d+"
  tiff_patterns: ["*.tif"]
  substitutions:
    - pattern: "exp_(\\d+)"
      repl: "experiment\\\\1"
```

### Animal ID-based Datasets
```yaml
dataset_discovery:
  pattern: "mouse_[A-Z]\\d{3}"  # Match mouse IDs like mouse_A123
  tiff_patterns: ["*.tif"]
  substitutions:
    - pattern: "mouse_([A-Z]\\d{3})"
      repl: "subject_\\\\1"
```

### Session-based Datasets
```yaml
dataset_discovery:
  pattern: "session_\\d{3}"  # Match session_001, session_002, etc.
  tiff_patterns: ["*.tif"]
  substitutions:
    - pattern: "session_(\\d{3})"
      repl: "s\\\\1"  # Convert to shorter format like s001
```

### Multi-level Directory Structure
```yaml
dataset_discovery:
  pattern: "subject_\\d+/session_\\d+"  # Match subject_1/session_1, etc.
  tiff_patterns: ["*.tif"]
  substitutions:
    - pattern: "subject_(\\d+)/session_(\\d+)"
      repl: "s\\\\1_s\\\\2"  # Convert to s1_s1 format
  exclude_patterns:
    - ".*/test/.*"  # Exclude test directories
    - ".*/backup/.*"  # Exclude backup directories
```

### Complex Pattern Matching
```yaml
dataset_discovery:
  pattern: "^(?:raw|processed)_\\d{8}_[A-Z]{2}"  # Match raw_20240315_AB or processed_20240315_AB
  tiff_patterns: ["*.tif"]
  substitutions:
    - pattern: "^(raw|processed)_(\\d{8})_([A-Z]{2})"
      repl: "\\\\2_\\\\3_\\\\1"  # Reorder to 20240315_AB_raw
```

The discovered datasets are automatically used in the Snakemake workflow.
