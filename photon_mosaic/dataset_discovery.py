"""
Dataset discovery module.

This module provides functions to discover datasets using regex patterns.
All filtering and transformations are handled through regex substitutions.
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union


def discover_datasets(
    base_path: Union[str, Path],
    pattern: str = ".*",
    exclude_patterns: Optional[List[str]] = None,
    substitutions: Optional[List[Dict[str, str]]] = None,
    tiff_patterns: list = ["*.tif"],
) -> Tuple[List[str], List[str], Dict[str, Dict[int, List[str]]], List[str]]:
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
    tiff_patterns : list, optional
        Glob pattern for TIFF files (default: '*.tif').

    Returns
    -------
    Tuple[List[str], List[str], Dict[str, Dict[int, List[str]]], List[str]]
        - List of original dataset names
        - List of transformed dataset names
        - Dictionary mapping original dataset names to their TIFF files by session
        - List of all TIFF files found

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

    datasets = sorted(datasets)
    original_datasets = sorted(original_datasets)

    # Discover TIFF files for each dataset
    tiff_files = {}
    tiff_files_flat = []

    for dataset in original_datasets:
        dataset_path = base_path_obj / dataset

        # check if there is at least one tiff in the dataset
        if not any(dataset_path.rglob("*.tif")):
            logging.info(f"No tiff files found in {dataset_path}")
            idx = datasets.index(dataset)
            datasets.pop(idx)
            original_datasets.pop(idx)
            continue

        # Initialize the dataset entry with all sessions
        tiff_files[dataset] = {}

        for session, tiff_pattern in enumerate(tiff_patterns):
            logging.debug(
                f"Searching for tiff files in {dataset_path} with pattern {tiff_pattern}"
            )
            files_found = sorted(
                [
                    f.name
                    for f in dataset_path.rglob(tiff_pattern)
                    if f.is_file()
                ]
            )

            if not files_found:
                logging.info(
                    f"No files found for pattern {tiff_pattern} in {dataset_path}"
                )
                # Initialize empty list for this session
                tiff_files[dataset][session] = []
            else:
                tiff_files[dataset][session] = files_found
                tiff_files_flat.extend(files_found)

    return original_datasets, datasets, tiff_files, tiff_files_flat
