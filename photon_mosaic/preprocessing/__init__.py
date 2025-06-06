"""
Preprocessing registry for photon-mosaic.

This module provides a registry for preprocessing steps that can be applied to image data.
Each preprocessing step is a function that takes an image array and returns a processed image array.
"""

from . import contrast
from . import derotation
from . import noop
from .registry import register_step, get_step, list_steps

__all__ = ["register_step", "get_step", "list_steps", "contrast", "derotation", "noop"]
