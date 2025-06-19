from pathlib import Path
from photon_mosaic.rules.preprocessing import run_preprocessing
from photon_mosaic.pathing import cross_platform_path
import re
import logging
import os

# Preprocessing rule
rule preprocessing:
    input:
        img=lambda wildcards: cross_platform_path(raw_data_base / datasets_old_names[int(wildcards.sub_idx)])
    output:
        processed=cross_platform_path(
            Path(processed_data_base).resolve()
            / "sub-{sub_idx}_{dataset}"
            / "ses-{ses_idx}"
            / "funcimg"
            / (f"{output_pattern}"+ "{tiff}")
        )
    params:
        dataset_folder=lambda wildcards: cross_platform_path(raw_data_base / datasets_old_names[int(wildcards.sub_idx)]),
        output_folder=lambda wildcards: cross_platform_path(
            Path(processed_data_base).resolve()
            / f"sub-{wildcards.sub_idx}_{datasets_new_names[int(wildcards.sub_idx)]}"
            / f"ses-{wildcards.ses_idx}"
            / "funcimg"
        ),
    wildcard_constraints:
        tiff="|".join(sorted(tiff_files_flat)),
        dataset="|".join(datasets_new_names),
    resources:
        **(slurm_config if config.get("use_slurm") else {}),
    run:
        from photon_mosaic.rules.preprocessing import run_preprocessing
        run_preprocessing(
            Path(params.output_folder),
            config["preprocessing"],
            Path(params.dataset_folder),
            ses_idx=int(wildcards.ses_idx),
            tiff_name=wildcards.tiff,
        )
