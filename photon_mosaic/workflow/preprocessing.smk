from pathlib import Path
from photon_mosaic.rules.preprocessing import run_preprocessing
import re


tiff_regex = "|".join(re.escape(name) for name in output_patterns)

# Preprocessing rule
rule preprocessing:
    input:
        tiffs=lambda wildcards: get_input_files(
            raw_data_base / datasets_old_names[int(wildcards.sub_idx)],
            config
        )
    output:
        processed="{processed_data}/sub-{sub_idx}_{dataset}/ses-{ses_idx}/funcimg/{tiff}"
    params:
        dataset_folder=lambda wildcards: str(raw_data_base / datasets_old_names[int(wildcards.sub_idx)])
    wildcard_constraints:
        tiff=tiff_regex
    run:
        from photon_mosaic.rules.preprocessing import run_preprocessing
        run_preprocessing(
            output.processed,
            config["preprocessing"],
            Path(params.dataset_folder),
        )
