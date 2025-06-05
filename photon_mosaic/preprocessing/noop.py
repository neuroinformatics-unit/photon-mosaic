"""
No-operation preprocessing step for photon-mosaic.

This module provides a function that returns the input data unchanged.
This is useful when preprocessing should be skipped.
"""

import numpy as np

from .registry import register_step


@register_step("noop")
def run(data: np.ndarray, **kwargs):
    """
    No-operation preprocessing step.

    Parameters
    ----------
    data : np.ndarray
        The image data.
    **kwargs : dict
        Additional arguments (ignored).
    """
    pass
