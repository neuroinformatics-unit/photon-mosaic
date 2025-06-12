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

Here's a detailed explanation of the configuration file structure:

```yaml
# ============================================================================
# Dataset Discovery Configuration
# ============================================================================
dataset_discovery:
  # Pattern to match dataset directories (e.g., date-based folders)
  pattern: "^2.*"
  
  # Patterns to exclude from processing
  exclude_patterns:
    - ".*_test$"    # Exclude test datasets
    - ".*_backup$"  # Exclude backup datasets
  
  # Transform dataset names using regex substitutions
  substitutions:
    # Remove underscores from names
    - pattern: "_"
      repl: ""
    # Take the first part of the path
    - pattern: "([^/]+)/.*"
      repl: '\\1'

# ============================================================================
# Data Paths Configuration
# ============================================================================
# Base directories for raw and processed data
raw_data_base: "/path/to/raw/"      # Directory containing raw imaging data
processed_data_base: "/path/to/processed/"  # Directory for processed outputs

# ============================================================================
# Suite2p Configuration
# ============================================================================
# Parameters for Suite2p processing pipeline
suite2p_ops:
  # Version and basic settings
  look_one_level_down: false
  fast_disk: []
  delete_bin: false
  mesoscan: false
  bruker: false
  bruker_bidirectional: false
  
  # File handling settings
  h5py: []
  h5py_key: 'data'
  nwb_file: ''
  nwb_driver: ''
  nwb_series: ''
  save_path0: ''
  save_folder: []
  subfolders: []
  move_bin: false
  save_mat: false
  save_NWB: false
  
  # Acquisition parameters
  nplanes: 1
  nchannels: 1
  functional_chan: 1
  tau: 1.0
  fs: 10.0
  force_sktiff: false
  frames_include: -1
  multiplane_parallel: false
  ignore_flyback: []
  
  # Registration settings
  do_registration: true
  two_step_registration: false
  keep_movie_raw: false
  nimg_init: 300
  batch_size: 500
  maxregshift: 0.1
  align_by_chan: 1
  reg_tif: false
  reg_tif_chan2: false
  subpixel: 10
  smooth_sigma_time: 0
  smooth_sigma: 1.15
  th_badframes: 1.0
  norm_frames: true
  force_refImg: false
  pad_fft: false
  nonrigid: true
  block_size: [128, 128]
  snr_thresh: 1.2
  maxregshiftNR: 5
  1Preg: false
  spatial_hp_reg: 42
  pre_smooth: 0
  spatial_taper: 40
  refImg_min_percentile: 1  # Minimum percentile for reference image selection (custom fork feature)
  
  # ROI detection settings
  roidetect: true
  spikedetect: true
  sparse_mode: true
  spatial_scale: 0
  connected: true
  nbinned: 5000
  max_iterations: 20
  threshold_scaling: 1.0
  max_overlap: 0.75
  high_pass: 100
  spatial_hp_detect: 25
  denoise: false
  anatomical_only: 0
  diameter: 0
  cellprob_threshold: 0.0
  flow_threshold: 1.5
  spatial_hp_cp: 0
  pretrained_model: 'cyto'
  
  # Neuropil extraction settings
  soma_crop: true
  neuropil_extract: true
  inner_neuropil_radius: 2
  min_neuropil_pixels: 350
  lam_percentile: 50.0
  allow_overlap: false
  
  # Classification settings
  use_builtin_classifier: false
  classifier_path: ''
  chan2_thres: 0.65
  
  # Baseline settings
  baseline: 'maximin'
  win_baseline: 60.0
  sig_baseline: 10.0
  prctile_baseline: 8.0
  neucoeff: 0.7

# ============================================================================
# Preprocessing Configuration
# ============================================================================
# Controls the preprocessing steps applied to the data
preprocessing:
  # Output file naming pattern for preprocessed files
  # Use {tiff_name} as a placeholder for the original tiff filename
  output_pattern: "preprocessed_{tiff_name}.tif"
  
  # Define preprocessing steps and their parameters
  steps:
    # Step 1: Derotation
    - name: derotation
      kwargs:
        # File pattern matching for input files
        glob_naming_pattern_tif: "*.tif"
        glob_naming_pattern_bin: "*.bin"
        path_to_stimulus_randperm: "/path/to/stimulus_randperm.npy"
    
    # Step 2: Contrast Enhancement
    - name: contrast
      kwargs:
        # File pattern matching for input files
        glob_naming_pattern_tif: "*.tif"
        # Contrast enhancement parameters
        percentile_low: 1    # Lower percentile for contrast stretching
        percentile_high: 99  # Upper percentile for contrast stretching

# ============================================================================
# SLURM Configuration
# ============================================================================
# Settings for SLURM job scheduling system
use_slurm: true
slurm:
  slurm_partition: "gpu"  # GPU partition for processing
  mem_mb: 32000          # Memory allocation in MB
  tasks: 1               # Number of tasks
  nodes: 1               # Number of nodes to use
```

## Key Parameters Explained

### Dataset Discovery
- `pattern`: Regular expression to identify dataset directories
- `exclude_patterns`: Patterns to skip during processing
- `substitutions`: Rules to transform dataset names

### Data Paths
- `raw_data_base`: Root directory for raw imaging data
- `processed_data_base`: Directory for processed outputs

### Suite2p Parameters
- **Acquisition**: Basic imaging parameters like sampling rate and number of planes
- **Registration**: Motion correction settings
- **ROI Detection**: Cell detection and segmentation parameters
- **Neuropil**: Neuropil extraction and decontamination settings
- **Classification**: Cell classification parameters
- **Baseline**: Baseline correction settings

### Preprocessing
- `output_pattern`: Template for output filenames
- `steps`: List of preprocessing steps with their parameters

### SLURM
- `use_slurm`: Enable/disable SLURM job scheduling
- `slurm_partition`: Compute partition to use
- `mem_mb`: Memory allocation per job
- `tasks`: Number of parallel tasks
- `nodes`: Number of compute nodes

## Tips for Configuration

1. **Start with Defaults**: Begin with the default configuration and modify only what you need
2. **Test with Dry Run**: Use `--dry-run` to verify your configuration
3. **Check Logs**: Review the generated logs in `derivatives/photon-mosaic/logs/`
4. **Version Control**: Keep your project-specific configs in version control
5. **Document Changes**: Note any configuration changes in your lab notebook or documentation
