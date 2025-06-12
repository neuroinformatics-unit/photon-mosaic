(user_guide/configuration)=
# Configuration

The configuration system in `photon-mosaic` is designed to be flexible and user-friendly. It allows you to customize the behavior of the pipeline at different levels.

## Configuration Files

### User Configuration
On first run, photon-mosaic will create a user config at `~/.photon_mosaic/config.yaml` if it does not exist. This serves as your default configuration.

### Project Configuration
You can also create a project-specific configuration file to override the defaults for a particular analysis.

## Configuration Sources

You can specify your configuration in three ways:

1. **User-wide defaults**: Edit the file at `~/.photon_mosaic/config.yaml` directly
2. **Command line overrides**: Override specific paths at runtime:
   ```bash
   photon-mosaic --raw_data_base /my/data --processed_data_base /my/processed --jobs 5
   ```
3. **Project-specific config**: Use a custom config file:
   ```bash
   photon-mosaic --config ./my/path/to/config.yaml --jobs 5
   ```

## Configuration Tracking

For reproducibility and debugging:
- The config used for each run (with any overrides) is exported to `derivatives/photon-mosaic/configs/YYYYMMDD_HHMMSS_config.yaml`
- Snakemake logs are saved to `derivatives/photon-mosaic/logs/YYYYMMDD_HHMMSS_snakemake.log`
- Both use timestamps (format: YYYYMMDD_HHMMSS) for easy tracking of different runs

## Configuration Structure

The configuration file is organized into several main sections. Here a simplified example:

```yaml
# Data paths
raw_data_base: "/path/to/raw/"
processed_data_base: "/path/to/processed/"

# Dataset discovery settings
dataset_discovery:
  pattern: "^2.*"  # Pattern to match dataset directories
  exclude_patterns: []  # Patterns to skip
  substitutions: []  # Rules to transform names

# Suite2p settings
suite2p_ops:
  # Some acquisition parameters
  nplanes: 1
  nchannels: 1
  fs: 10.0
  tau: 1.0
  # any other Suite2p parameters can be added here

  # Cell detection
  anatomical_only: 0  # Set > 0 to use Cellpose

# Preprocessing steps
preprocessing:
  steps:
    - name: contrast
      kwargs:
        glob_naming_pattern_tif: "*.tif"
        percentile_low: 1
        percentile_high: 99

# SLURM settings
use_slurm: true
slurm:
  partition: "gpu"
  mem_mb: 32000
  tasks: 1
  nodes: 1
```
For the full configuration file, see [photon_mosaic/workflow/config.yaml](https://github.com/neuroinformatics-unit/photon-mosaic/blob/main/photon_mosaic/workflow/config.yaml) or the YAML file in `~/.photon-mosaic/config.yaml` generated on first run.

## Key Parameters Explained

### Dataset Discovery
- `pattern`: Regular expression to identify dataset directories
- `exclude_patterns`: Patterns to skip during processing
- `substitutions`: Rules to transform dataset names

### Data Paths
- `raw_data_base`: Root directory for raw imaging data
- `processed_data_base`: Directory for processed outputs

### Suite2p Parameters
For a complete list of all available Suite2p parameters and their descriptions, please refer to the [official Suite2p documentation](https://suite2p.readthedocs.io/en/latest/settings.html).

#### Registration Settings
Our custom fork of Suite2p includes additional parameters for registration:
- `refImg_min_percentile`: Minimum percentile for reference image selection (default: 1). Controls the lower bound of contrast normalization. Higher values improve registration for low SNR datasets like three-photon imaging.
- `refImg_max_percentile`: Maximum percentile for reference image selection (default: 99). Controls the upper bound of contrast normalization. Together with `refImg_min_percentile`, defines the contrast range for registration.

#### Cell Detection
To use Cellpose for cell detection, set `anatomical_only` to a value greater than 0 in your configuration. For example:

```yaml
suite2p_ops:
  anatomical_only: 4  # Use maximum projection image for cell detection
```

### Preprocessing
For detailed information about preprocessing steps and their configuration, see the [preprocessing documentation](preprocessing.md).

### SLURM
- `use_slurm`: Enable/disable SLURM job scheduling
- `slurm_partition`: Compute partition to use
- `mem_mb`: Memory allocation per job
- `tasks`: Number of parallel tasks
- `nodes`: Number of compute nodes

For more details about SLURM configuration options, see the [Snakemake SLURM executor plugin documentation](https://github.com/snakemake/snakemake-executor-plugin-slurm).

## Tips for Configuration

1. **Start with Defaults**: Begin with the default configuration and modify only what you need
2. **Test with Dry Run**: Use `--dry-run` to verify your configuration
3. **Check Logs**: Review the generated logs in `derivatives/photon-mosaic/logs/`
4. **Version Control**: Keep your project-specific configs in version control
5. **Document Changes**: Note any configuration changes in your lab notebook or documentation
