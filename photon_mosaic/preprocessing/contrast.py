"""
Contrast enhancement module.

This module provides functions to enhance the contrast of image data.
"""

from skimage.exposure import equalize_adapthist
from tifffile import imwrite

from photon_mosaic.preprocessing.registry import register_step


@register_step("contrast")
def run(data, clip_limit=0.01, output_path=None, **kwargs):
    """
    Enhance the contrast of image data using CLAHE.

    Parameters
    ----------
    data : numpy.ndarray
        Input image data.
    clip_limit : float, optional
        Clipping limit for CLAHE, by default 0.01.
    output_path : str, optional
        Path to save the enhanced image. If None, the image is not saved.
    **kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    numpy.ndarray
        Enhanced image data.
    """
    # Normalize data to [0, 1] range
    data_norm = (data - data.min()) / (data.max() - data.min())

    # Apply CLAHE
    enhanced = equalize_adapthist(data_norm, clip_limit=clip_limit)

    # Scale back to original range
    enhanced = enhanced * (data.max() - data.min()) + data.min()

    # Save the enhanced image if output path is provided
    if output_path is not None:
        imwrite(output_path, enhanced.astype(data.dtype))
