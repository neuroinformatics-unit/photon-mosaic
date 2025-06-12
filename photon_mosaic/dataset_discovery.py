"""
Dataset discovery module.

This module provides functions to discover datasets using regex patterns.
All filtering and transformations are handled through regex substitutions.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Union


def discover_datasets(
    base_path: Union[str, Path],
    pattern: str = ".*",
    exclude_patterns: Optional[List[str]] = None,
    substitutions: Optional[List[Dict[str, str]]] = None,
) -> List[str]:
    """
    Discover datasets in a directory using regex patterns.

    Parameters
    ----------
    base_path : str or Path
        Base path to search for datasets.
    pattern : str, optional
        Regex pattern to match dataset names, defaults to ".*"
        (all directories).
    exclude_patterns : List[str], optional
        List of regex patterns for datasets to exclude.
    substitutions : List[Dict[str, str]], optional
        List of regex substitution pairs to transform dataset names.
        Each dict should have 'pattern' and 'repl' keys for re.sub().

    Returns
    -------
    List[str]
        List of discovered dataset names.

    Examples
    --------
    >>> # Find all datasets starting with '2' (e.g., date-based)
    >>> datasets = discover_datasets("/path/to/data", pattern="^2.*")

    >>> # Exclude test datasets and remove underscores
    >>> datasets = discover_datasets(
    ...     "/path/to/data",
    ...     pattern="^2.*",
    ...     exclude_patterns=[".*_test$"],
    ...     substitutions=[{"pattern": "_", "repl": ""}],
    ... )

    >>> # Complex transformations using regex groups
    >>> datasets = discover_datasets(
    ...     "/path/to/data",
    ...     pattern=".*",
    ...     substitutions=[
    ...         # Convert "exp_001" to "experiment001"
    ...         {"pattern": "exp_(\d+)", "repl": r"experiment\1"},
    ...         # Remove any trailing _v1, _v2, etc.
    ...         {"pattern": "_v\d+$", "repl": ""},
    ...     ],
    ... )
    """
    # Convert base_path to Path if it's a string
    base_path_obj = (
        Path(base_path) if isinstance(base_path, str) else base_path
    )

    # Find all directories matching the pattern
    datasets = [
        d.name
        for d in base_path_obj.iterdir()
        if d.is_dir() and re.match(pattern, d.name)
    ]

    # Apply exclusion patterns
    if exclude_patterns:
        for exclude in exclude_patterns:
            datasets = [ds for ds in datasets if not re.match(exclude, ds)]

    # Apply regex substitutions
    if substitutions:
        for sub in substitutions:
            datasets = [
                re.sub(sub["pattern"], sub["repl"], ds) for ds in datasets
            ]

    return sorted(datasets)
