"""
Preprocessing rule for photon-mosaic Snakemake workflow.

This file contains the rule for preprocessing image data before running Suite2P.
"""

from pathlib import Path
from photon_mosaic.rules.preprocessing import run_preprocessing

rule preprocessing:
    input:
        tiff=lambda wildcards: str(tiff_paths[wildcards.datasets][int(wildcards.tiff_index)])
    output:
        processed=f"{processed_data_base}/sub-{{index}}_{{datasets}}/ses-0/funcimg/{{tiff_name}}.tif"
    params:
        dataset_folder=lambda wildcards: raw_data_base / wildcards.datasets,
    resources:
        **(slurm_config if slurm_config.get("use_slurm") else {}),
    run:
        from photon_mosaic.rules.preprocessing_run import run_preprocessing
        run_preprocessing(
            input["tiff"],
            output["processed"],
            Path(params["dataset_folder"]),
            config["preprocessing_ops"],
        )

# Expand the preprocessing rule for all tiff files
expand(
    "{processed_data_base}/sub-{index}_{datasets}/ses-0/funcimg/{tiff_name}.tif",
    processed_data_base=processed_data_base,
    zip,
    index=[i for i, _ in enumerate(datasets)],
    datasets=datasets,
    tiff_name=[tiff["tiff_name"] for tiff in tiff_files],
)
