"""
Dataset discovery module.

This module provides functions to discover datasets using regex patterns.
All filtering and transformations are handled through regex substitutions.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union


def discover_datasets(
    base_path: Union[str, Path],
    pattern: str = ".*",
    exclude_patterns: Optional[List[str]] = None,
    substitutions: Optional[List[Dict[str, str]]] = None,
    tiff_pattern: str = "*.tif",
) -> Tuple[List[str], List[str], Dict[str, List[str]]]:
    """
    Discover datasets and their TIFF files in a directory using regex patterns.

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
    tiff_pattern : str, optional
        Glob pattern for TIFF files (default: '*.tif').

    Returns
    -------
    Tuple[List[str], List[str], Dict[str, List[str]]]
        - List of original dataset names
        - List of transformed dataset names
        - Dictionary mapping original dataset names to their TIFF files

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

    # Store original dataset names
    original_datasets = datasets.copy()

    # Apply regex substitutions to get new names
    if substitutions:
        for sub in substitutions:
            datasets = [
                re.sub(sub["pattern"], sub["repl"], ds) for ds in datasets
            ]

    # Discover TIFF files for each dataset
    tiff_files = {}
    for dataset in original_datasets:
        dataset_path = base_path_obj / dataset
        tiff_files[dataset] = sorted([
            f.name for f in dataset_path.glob(tiff_pattern)
        ])

    return sorted(original_datasets), sorted(datasets), tiff_files

