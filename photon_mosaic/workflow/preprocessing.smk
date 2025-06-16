from pathlib import Path
from photon_mosaic.rules.preprocessing import run_preprocessing
import re


#  Necessary for the suite2p rule to understand the input mapping
tiff_regex = "|".join(re.escape(name) for name in output_patterns)

# Preprocessing rule
rule preprocessing:
    input:
        tiff=lambda wildcards: get_input_files(
            dataset_folder = raw_data_base / datasets_old_names[int(wildcards.sub_idx)],
            config = config,
            ses_idx = int(wildcards.ses_idx),
        )
    output:
        processed=str(
            Path("{processed_data}")
            / "sub-{sub_idx}_{dataset}"
            / "ses-{ses_idx}"
            / "funcimg"
            / "{tiff}"
        )
    params:
        dataset_folder=lambda wildcards: str(raw_data_base / datasets_old_names[int(wildcards.sub_idx)]),
        output_folder=lambda wildcards: str(
            Path(processed_data_base)
            / f"sub-{wildcards.sub_idx}_{datasets_new_names[int(wildcards.sub_idx)]}"
            / f"ses-{wildcards.ses_idx}"
            / "funcimg"
        ),
    wildcard_constraints:
        tiff=tiff_regex
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
