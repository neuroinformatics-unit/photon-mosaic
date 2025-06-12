"""
Registry for preprocessing steps.

This module provides the registry functionality for preprocessing steps.
"""

from typing import Callable, Dict, Optional

# Registry to store preprocessing steps
_PREPROCESSING_STEPS: Dict[str, Callable] = {}


def register_step(name: Optional[str] = None) -> Callable:
    """
    Decorator to register a preprocessing step.

    Parameters
    ----------
    name : str, optional
        The name to register the step under. If not provided, the function name
        is used.

    Returns
    -------
    Callable
        The decorated function.
    """

    def decorator(func: Callable) -> Callable:
        step_name = name or func.__name__
        _PREPROCESSING_STEPS[step_name] = func
        return func

    return decorator


def get_step(name: str) -> Callable:
    """
    Get a preprocessing step by name.

    Parameters
    ----------
    name : str
        The name of the preprocessing step.

    Returns
    -------
    Callable
        The preprocessing step function.

    Raises
    ------
    KeyError
        If the step is not found in the registry.
    """
    if name not in _PREPROCESSING_STEPS:
        raise KeyError(f"Preprocessing step '{name}' not found in registry")
    return _PREPROCESSING_STEPS[name]


def list_steps() -> Dict[str, Callable]:
    """
    List all registered preprocessing steps.

    Returns
    -------
    Dict[str, Callable]
        A dictionary mapping step names to their functions.
    """
    return _PREPROCESSING_STEPS.copy()
