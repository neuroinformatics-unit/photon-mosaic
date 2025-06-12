# Preprocessing

The preprocessing module in photon-mosaic provides a flexible system for applying preprocessing steps to image data. Each preprocessing step is a function that takes an image array and returns a processed image array.

## Available Preprocessing Steps

### Contrast Enhancement

The contrast enhancement step uses percentile-based contrast stretching to improve image contrast:

```python
from photon_mosaic.preprocessing import get_step

# Get the contrast enhancement step
contrast = get_step("contrast")

# Apply contrast enhancement
enhanced = contrast(
    dataset_folder="path/to/dataset",
    output_folder="path/to/output",
    ses_idx=0,
    glob_naming_pattern_tif="*.tif",
    percentile_low=1,
    percentile_high=99
)
```

The contrast enhancement works by:
1. Finding the pixel values at the specified percentiles
2. Stretching the image intensity to use the full range between these values
3. Saving the enhanced image with the prefix "enhanced_" in the output folder

**Note**: If you are using the `refImg_min_percentile` option in Suite2p (available in our custom fork), you may not need to perform contrast enhancement as a separate preprocessing step. This option automatically handles contrast normalization during the registration process. See the Suite2p configuration section for more details.

### Derotation

The derotation step handles image derotation using the derotation package:

```python
from photon_mosaic.preprocessing import get_step

# Get the derotation step
derotation = get_step("derotation")

# Apply derotation
derotated = derotation(
    dataset_folder="path/to/dataset",
    output_folder="path/to/output",
    ses_idx=0,
    glob_naming_pattern_tif="*.tif",
    glob_naming_pattern_bin="*.bin"
)
```

### No-Operation (Noop)

The noop step is useful when you want to skip preprocessing for certain files. It either returns the input data unchanged or copies the input file to the output directory:

```python
from photon_mosaic.preprocessing import get_step

# Get the noop step
noop = get_step("noop")

# Option 1: Pass data directly
processed = noop(data=image_array)

# Option 2: Copy file without modification
noop(
    dataset_folder="path/to/dataset",
    output_folder="path/to/output",
    ses_idx=0,
    glob_naming_pattern_tif="*.tif"
)
```

For detailed parameter documentation, please refer to the [API Reference](api_reference.html).

## Using Preprocessing in Configuration

In your configuration file (e.g., `config.yaml`), you can specify which preprocessing steps to apply and their parameters:

```yaml
preprocessing:
  steps:
    - name: derotation
      kwargs:
        glob_naming_pattern_tif: "*.tif"
        glob_naming_pattern_bin: "*.bin"
        path_to_stimulus_randperm: "/path/to/stimulus_randperm.npy"
    - name: contrast
      kwargs:
        glob_naming_pattern_tif: "*.tif"
        percentile_low: 1
        percentile_high: 99
    - name: noop
      kwargs:
        glob_naming_pattern_tif: "*.tif"
```

## Adding New Preprocessing Steps

To add a new preprocessing step:

1. Create a new module in the `photon_mosaic/preprocessing` directory
2. Define a function that takes a numpy array as input and returns a processed numpy array
3. Use the `@register_step` decorator to register the function
4. Import the module in `__init__.py` to ensure it's registered

Example:
```python
from photon_mosaic.preprocessing import register_step

@register_step("my_step")
def run(data, **kwargs):
    # Process the data
    return processed_data
```

## Best Practices

- Each preprocessing step should handle both 2D (single image) and 3D (stack of images) inputs
- Preprocessing steps should preserve the data type of the input
- Document the parameters and return values of each preprocessing step
- Include tests for each preprocessing step
