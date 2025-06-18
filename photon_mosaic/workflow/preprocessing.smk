from pathlib import Path
from photon_mosaic.rules.preprocessing import run_preprocessing
import re
import logging

tiff_regex = "|".join(TIFF_NAMES)
logging.debug(f"TIFF regex: {tiff_regex}")

# Preprocessing rule
rule preprocessing:
    input:
        img=lambda wildcards: get_input_files(
            dataset_folder = raw_data_base / datasets_old_names[int(wildcards.sub_idx)],
            config = config,
            ses_idx = int(wildcards.ses_idx),
        )
    output:
        processed=str(
            Path(processed_data_base).resolve()
            / "sub-{sub_idx}_{dataset}"
            / "ses-{ses_idx}"
            / "funcimg"
            / "{tiff}"
        )
    params:
        dataset_folder=lambda wildcards: str(raw_data_base / datasets_old_names[int(wildcards.sub_idx)]),
        output_folder=lambda wildcards: str(
            Path(processed_data_base).resolve()
            / f"sub-{wildcards.sub_idx}_{datasets_new_names[int(wildcards.sub_idx)]}"
            / f"ses-{wildcards.ses_idx}"
            / "funcimg"
        ),
    wildcard_constraints:
        tiff=tiff_regex,
        dataset="|".join(datasets_new_names)
    resources:
        **(slurm_config if config.get("use_slurm") else {}),
    run:
        from photon_mosaic.rules.preprocessing import run_preprocessing
        run_preprocessing(
            Path(params.output_folder),
            config["preprocessing"],
            Path(params.dataset_folder),
            ses_idx=int(wildcards.ses_idx),
        )
