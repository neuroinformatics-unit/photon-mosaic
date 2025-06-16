rule suite2p:
    input:
        tiff=lambda wildcards: str(
            Path(processed_data_base).resolve()
            / f"sub-{wildcards.sub_idx}_{datasets_new_names[int(wildcards.sub_idx)]}"
            / f"ses-{wildcards.ses_idx}"
            / "funcimg"
            / f"{output_patterns[int(wildcards.ses_idx)]}"
        )
    output:
        F=str(
            Path(processed_data_base).resolve()
            / "sub-{sub_idx}_{dataset}"
            / "ses-{ses_idx}"
            / "funcimg"
            / "suite2p"
            / "plane0"
            / "F.npy"
        ),
        bin=str(
            Path(processed_data_base).resolve()
            / "sub-{sub_idx}_{dataset}"
            / "ses-{ses_idx}"
            / "funcimg"
            / "suite2p"
            / "plane0"
            / "data.bin"
        )
    params:
        dataset_folder=lambda wildcards: str(
            Path(processed_data_base).resolve()
            / f"sub-{wildcards.sub_idx}_{datasets_new_names[int(wildcards.sub_idx)]}"
            / f"ses-{wildcards.ses_idx}"
            / "funcimg"
        )
    wildcard_constraints:
        dataset="|".join(datasets_new_names)
    resources:
        **(slurm_config if config.get("use_slurm") else {}),
    run:
        from photon_mosaic.rules.suite2p_run import run_suite2p
        from pathlib import Path

        # Ensure all paths are properly resolved
        input_path = Path(input.tiff).resolve()
        output_path = Path(output.F).resolve()
        dataset_folder = Path(params.dataset_folder).resolve()

        run_suite2p(
            str(output_path),
            dataset_folder,
            config["suite2p_ops"],
        )
