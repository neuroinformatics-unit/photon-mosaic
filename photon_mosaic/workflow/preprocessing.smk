"""
Preprocessing rule for photon-mosaic Snakemake workflow.

This file contains the rule for preprocessing image data before running Suite2P.
"""

from photon_mosaic.rules.preprocessing import run_preprocessing, get_output_pattern

rule preprocessing:
    input:
        tiff=lambda wildcards: tiff_paths[wildcards.datasets],
    output:
        tiff=lambda wildcards: f"{processed_data_base}/sub-{wildcards.index}_{wildcards.datasets}/ses-0/funcimg/{get_output_pattern(wildcards.tiff_name, config, tiff_paths)}",
    params:
        dataset_folder=lambda wildcards: raw_data_base / wildcards.datasets,
        tiff_name=lambda wildcards: Path(wildcards.tiff).stem,
    resources:
        **(slurm_config if slurm_config.get("use_slurm") else {}),
    run:
        run_preprocessing(
            input.tiff,
            output.tiff,
            config,
            dataset_folder=str(params.dataset_folder)
        )
