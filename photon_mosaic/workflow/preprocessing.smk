"""
Preprocessing Module

This Snakefile module handles the preprocessing step of the photon mosaic pipeline.
It processes raw TIFF files from discovered datasets and applies preprocessing
operations defined in the configuration.

The preprocessing rule:
- Takes raw TIFF files as input from the dataset discovery
- Applies preprocessing operations (defined in config["preprocessing"])
- Outputs processed files in a standardized NeuroBlueprint format
- Supports SLURM cluster execution with configurable resources

Input: Raw TIFF files from discovered datasets
Output: Preprocessed TIFF files organized by subject/session
"""

from pathlib import Path
from photon_mosaic.rules.preprocessing import run_preprocessing
from photon_mosaic.pathing import cross_platform_path
import re
import logging
import os

# Preprocessing rule
rule preprocessing:
    input:
        img=lambda wildcards: cross_platform_path(raw_data_base / discoverer.original_datasets[discoverer.transformed_datasets.index(wildcards.subject_name)])
    output:
        processed=cross_platform_path(
            Path(processed_data_base).resolve()
            / "{subject_name}"
            / "{session_name}"
            / "funcimg"
            / (f"{output_pattern}"+ "{tiff}")
        )
    params:
        dataset_folder=lambda wildcards: cross_platform_path(raw_data_base / discoverer.original_datasets[discoverer.transformed_datasets.index(wildcards.subject_name)]),
        output_folder=lambda wildcards: cross_platform_path(
            Path(processed_data_base).resolve()
            / wildcards.subject_name
            / wildcards.session_name
            / "funcimg"
        ),
        ses_idx=lambda wildcards: int(wildcards.session_name.split("_")[0].replace("ses-", "")),
    wildcard_constraints:
        tiff="|".join(sorted(discoverer.tiff_files_flat)) if discoverer.tiff_files_flat else "dummy",
        subject_name="|".join(discoverer.transformed_datasets),
        session_name="|".join([discoverer.get_session_name(i, session_idx) for i in range(len(discoverer.transformed_datasets))
                              for session_idx in discoverer.tiff_files[discoverer.original_datasets[i]].keys()]),
    resources:
        **(slurm_config if config.get("use_slurm") else {}),
    run:
        from photon_mosaic.rules.preprocessing import run_preprocessing
        run_preprocessing(
            Path(params.output_folder),
            config["preprocessing"],
            Path(params.dataset_folder),
            ses_idx=int(params.ses_idx),
            tiff_name=wildcards.tiff,
        )
