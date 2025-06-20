"""
Suite2p Analysis Module

This Snakefile module handles the Suite2p analysis step of the photon mosaic pipeline.
It takes preprocessed TIFF files and runs Suite2p to extract neural activity traces.

The suite2p rule:
- Takes preprocessed TIFF files as input (from preprocessing step)
- Runs Suite2p analysis with parameters from config["suite2p_ops"]
- Outputs F.npy (fluorescence traces) and data.bin (binary data) files
- Supports SLURM cluster execution with configurable resources

Input: Preprocessed TIFF files from the preprocessing step
Output: Suite2p analysis results (F.npy, data.bin) in suite2p/plane0/ directory
"""

import re
from photon_mosaic.pathing import cross_platform_path

rule suite2p:
    input:
        tiffs=lambda wildcards: [
            cross_platform_path(
                Path(processed_data_base).resolve()
                / f"sub-{wildcards.sub_idx}_{datasets_new_names[int(wildcards.sub_idx)]}"
                / f"ses-{wildcards.ses_idx}"
                / "funcimg"
                / f"{output_pattern}{tiff_name}"
            )
            for tiff_name in tiff_files_map[int(wildcards.sub_idx)][int(wildcards.ses_idx)]
        ],
    output:
        F=cross_platform_path(
            Path(processed_data_base).resolve()
            / "sub-{sub_idx}_{dataset}"
            / "ses-{ses_idx}"
            / "funcimg"
            / "suite2p"
            / "plane0"
            / "F.npy"
        ),
        bin=cross_platform_path(
            Path(processed_data_base).resolve()
            / "sub-{sub_idx}_{dataset}"
            / "ses-{ses_idx}"
            / "funcimg"
            / "suite2p"
            / "plane0"
            / "data.bin"
        )
    params:
        dataset_folder=lambda wildcards: cross_platform_path(
            Path(processed_data_base).resolve()
            / f"sub-{wildcards.sub_idx}_{datasets_new_names[int(wildcards.sub_idx)]}"
            / f"ses-{wildcards.ses_idx}"
            / "funcimg"
        ),
    wildcard_constraints:
        dataset="|".join(datasets_new_names),
    resources:
        **(slurm_config if config.get("use_slurm") else {}),
    run:
        from photon_mosaic.rules.suite2p_run import run_suite2p
        from pathlib import Path

        # Ensure all paths are properly resolved
        input_paths = [Path(tiff).resolve() for tiff in input.tiffs]
        output_path = Path(output.F).resolve()
        dataset_folder = Path(params.dataset_folder).resolve()

        run_suite2p(
            str(output_path),
            dataset_folder,
            config["suite2p_ops"],
        )
