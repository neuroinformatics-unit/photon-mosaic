# Preprocessing

The preprocessing module in photon-mosaic provides a flexible system for applying preprocessing steps to image data. Each preprocessing step is a function that takes an image array and returns a processed image array.

## Available Preprocessing Steps

### Contrast Enhancement

The contrast enhancement step uses CLAHE (Contrast Limited Adaptive Histogram Equalization) to improve image contrast:

```python
from photon_mosaic.preprocessing import get_step

# Get the contrast enhancement step
contrast = get_step("contrast_enhancement")

# Apply contrast enhancement
enhanced = contrast(data, clip_limit=2.0)
```

Parameters:
- `data`: Input image data (numpy array)
- `clip_limit`: Threshold for contrast limiting (default: 2.0)

### Derotation

The derotation step handles image derotation using the derotation package:

```python
from photon_mosaic.preprocessing import get_step

# Get the derotation step
derotation = get_step("derotation")

# Apply derotation
derotated = derotation(data, dataset_folder="path/to/dataset", output_folder="path/to/output")
```

Parameters:
- `data`: Input image data (numpy array)
- `dataset_folder`: Path to the dataset folder
- `output_folder`: Path to the output folder

## Using Preprocessing in Configuration

In your configuration file (e.g., `config.yaml`), you can specify which preprocessing steps to apply and their parameters:

```yaml
preprocessing:
  steps:
    - name: derotation
      kwargs:
        dataset_folder: "path/to/dataset"
        output_folder: "path/to/output"
    - name: contrast_enhancement
      kwargs:
        clip_limit: 2.0
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
