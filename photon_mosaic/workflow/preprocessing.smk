"""
Preprocessing rule for photon-mosaic Snakemake workflow.

This file contains the rule for preprocessing image data before running Suite2P.
"""

from pathlib import Path
from photon_mosaic.rules.preprocessing import run_preprocessing, get_output_pattern

rule preprocessing:
    input:
        tiff=lambda wildcards: tiff_paths[wildcards.datasets],
    output:
        tiff=f"{processed_data_base}/sub-{{index}}_{{datasets}}/ses-0/funcimg/{{tiff_name}}.tif",
    params:
        dataset_folder=lambda wildcards: raw_data_base / wildcards.datasets,
        tiff_name=lambda wildcards: Path(wildcards.tiff).stem,
        output_pattern=lambda wildcards: get_output_pattern(Path(wildcards.tiff).stem, config, tiff_paths),
    resources:
        **(slurm_config if slurm_config.get("use_slurm") else {}),
    run:
        # Run preprocessing
        run_preprocessing(
            input.tiff,
            output.tiff,
            config,
            dataset_folder=str(params.dataset_folder)
        )
