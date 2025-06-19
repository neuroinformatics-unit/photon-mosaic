# Preprocessing

The preprocessing module in photon-mosaic provides a flexible system for applying preprocessing steps to image data. Each preprocessing step is a function that takes an image array and returns a processed image array.

## Available Preprocessing Steps

### Contrast Enhancement

The contrast enhancement step uses percentile-based contrast stretching to improve image contrast.

The contrast enhancement works by:
1. Finding the pixel values at the specified percentiles
2. Stretching the image intensity to use the full range between these values
3. Saving the enhanced image with the prefix "enhanced_" in the output folder

**Note**: If you are using the additional options in Suite2p (available in our custom fork), you may not need to perform contrast enhancement as a separate preprocessing step. This option automatically handles contrast normalization during the registration process. See the Suite2p configuration section for more details.

### Derotation

The derotation step handles image derotation using the derotation package. Please refer to the [derotation package documentation](https://derotation.neuroinformatics.dev/) for more details.

### No-Operation (Noop)

The noop step simply copies the input files to the output directory without any modification. Use this when you want to skip preprocessing but still need the files in the output directory.

## Configuration

To use preprocessing in your configuration, add a `preprocessing` section to your configuration file. Each step is defined as a dictionary with a `name` and optional `kwargs` for step-specific parameters.

### Basic Configuration

```yaml
preprocessing:
  output_pattern: "{tiff_name}.tif"
  steps:
    - name: contrast
      kwargs:
        glob_naming_pattern_tif: "*.tif"
        percentile_low: 1
        percentile_high: 99
```

### Step-Specific Parameters

Each preprocessing step can have its own configuration parameters:

```yaml
preprocessing:
  output_pattern: "{tiff_name}.tif"
  steps:
    - name: derotation
      kwargs:
        glob_naming_pattern_tif: "*.tif"
        glob_naming_pattern_bin: "*.bin"
        path_to_stimulus_randperm: "/path/to/stimulus_randperm.npy"
    - name: contrast
      kwargs:
        glob_naming_pattern_tif: "derotated_*.tif"
        percentile_low: 1
        percentile_high: 99
```

See the [API Reference](api_reference.html) for more details on the available parameters.

## Chaining and Selecting Steps

### Important Note on Execution

Preprocessing steps are executed sequentially in a single job, in the exact order they appear in your configuration file. Each step saves its output to disk and the next step must be configured to find these output files. Ideed, in the example above, the contrast step is configured to find the output files from the derotation step.

### File Naming Conventions

Each step has its own output naming convention:
- Derotation step: Uses the derotation package's naming convention, "derotated_<optional_suffix>.tif"
- Contrast step: Adds "enhanced_" prefix to input filenames
- Noop step: Copies files with original names

### Selecting Steps

To run only specific steps, comment out the ones you don't want or delete them.

```yaml
preprocessing:
  steps:
    # - name: derotation  # Commented out to skip
    #   kwargs:
    #     glob_naming_pattern_tif: "*.tif"
    - name: contrast
      kwargs:
        glob_naming_pattern_tif: "*.tif"
```

To skip all preprocessing:
```yaml
preprocessing:
  steps:
    - name: noop
      kwargs:
        glob_naming_pattern_tif: "*.tif"
```

## Adding New Preprocessing Steps

To add a new preprocessing step:

1. Create a new module in the `photon_mosaic/preprocessing` directory (e.g., `my_step.py`)
2. Define a function named `run` that takes the required parameters and processes the data
3. Import the module in `__init__.py` to make it available

Example:
```python
# photon_mosaic/preprocessing/my_step.py
def run(dataset_folder, output_folder, ses_idx, **kwargs):
    # Process the data
    # Save results to output_folder
    pass
```

The step will be automatically available by its module name (e.g., `my_step`) in your configuration.
