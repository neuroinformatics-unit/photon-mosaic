rule setup:
    input:
        datasets_path="/nfs/winstor/margrie/SimonWeiler/RawData/Invivo_imaging/3photon_rotation/shared/",
        writing_path="/ceph/margrie/laura/cimaut/",
    output: "setup_output.txt"
    shell: "python calcium_imaging_automation/core/rules/setup.py {input.datasets_path} {input.writing_path} --folder_read_pattern '2*' --file_read_pattern 'rotation_00001.tif' --file_read_pattern '*.bin' > {output}"
