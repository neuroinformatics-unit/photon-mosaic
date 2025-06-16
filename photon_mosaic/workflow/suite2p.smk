rule suite2p:
    input:
        tiff=lambda wildcards: str(
            Path(processed_data_base)
            / f"sub-{wildcards.sub_idx}_{datasets_new_names[int(wildcards.sub_idx)]}"
            / f"ses-{wildcards.ses_idx}"
            / "funcimg"
            / f"{output_patterns[int(wildcards.ses_idx)]}"
        )
    output:
        F=lambda wildcards, processed_data_base, dataset, sub_idx, ses_idx: str(
            Path(processed_data_base)
            / f"sub-{sub_idx}_{dataset}"
            / f"ses-{ses_idx}"
            / "funcimg"
            / "suite2p"
            / "plane0"
            / "F.npy"
        ),
        bin=lambda wildcards, processed_data_base, dataset, sub_idx, ses_idx: str(
            Path(processed_data_base)
            / f"sub-{sub_idx}_{dataset}"
            / f"ses-{ses_idx}"
            / "funcimg"
            / "suite2p"
            / "plane0"
            / "data.bin"
        )
    params:
        dataset_folder=lambda wildcards, processed_data_base, dataset, sub_idx, ses_idx: str(
            Path(processed_data_base)
            / f"sub-{sub_idx}_{dataset}"
            / f"ses-{ses_idx}"
            / "funcimg"
        )
    resources:
        **(slurm_config if config.get("use_slurm") else {}),
    run:
        from photon_mosaic.rules.suite2p_run import run_suite2p
        run_suite2p(
            output.F,
            Path(params.dataset_folder),
            config["suite2p_ops"],
        )
