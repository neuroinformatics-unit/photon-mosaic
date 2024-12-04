rule setup:
    input:
        datasets_path="/nfs/winstor/margrie/SimonWeiler/RawData/Invivo_imaging/3photon_rotation/shared/",
        writing_path="/ceph/margrie/laura/cimaut/",
    output: "setup_output.txt"
    shell: "python calcium_imaging_automation/core/rules/setup.py {input.datasets_path} {input.writing_path} --folder_read_pattern '2*' --file_read_pattern 'rotation_00001.tif' --file_read_pattern '*.bin' > {output}"

import pandas as pd

paths = pd.read_csv("datasets.csv")

rule all:
    input:
        expand("preprocess_output_{index}.txt", index=paths["index"])

rule preprocess:
    input:
        lambda wildcards: paths.loc[int(wildcards.index), "read_dataset_path"],
        lambda wildcards: paths.loc[int(wildcards.index), "write_dataset_path"],
    output:
        "preprocess_output_{index}.txt"
    params:
        index=lambda wildcards: wildcards.index
    script:
        "calcium_imaging_automation/core/rules/preprocess.py"
