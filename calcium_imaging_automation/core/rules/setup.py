import argparse
import shutil
from pathlib import Path

import pandas as pd

from calcium_imaging_automation.core.reader import ReadAquiredData
from calcium_imaging_automation.core.writer import DatashuttleWrapper
from snakemake.script import snakemake


try:
    read_dataset_path = Path(snakemake.input[0])
    write_dataset_path = Path(snakemake.input[1])
    folder_read_pattern = snakemake.params.folder_read_pattern
    file_read_pattern = snakemake.params.file_read_pattern
    
    output = snakemake.output[0]

    try:
        shutil.rmtree("/ceph/margrie/laura/cimaut/derivatives/")
        shutil.rmtree("/ceph/margrie/laura/cimaut/submitit/")
    except FileNotFoundError:
        print("No derivatives folder found")

    print(f"Reading data from {read_dataset_path}")

    reader = ReadAquiredData(
        read_dataset_path,
        folder_read_pattern,
        file_read_pattern,
    )
    print(f"Found {len(reader.datasets_paths)} datasets.")

    number_of_tiffs = reader.max_session_number(filetype="tif")
    print(f"Max of tiffs found: {number_of_tiffs}")

    writer = DatashuttleWrapper(write_dataset_path)
    writer.create_folders(reader.dataset_names, session_number=number_of_tiffs)
    print("Folders created")

    datasets = pd.DataFrame(
        {
            "read_dataset_path": reader.datasets_paths,
            "write_dataset_path": [
                writer.get_dataset_path(dt.stem)
                for dt in reader.datasets_paths
            ],
        }
    )
    datasets.to_csv(output, index=True, index_label="index")
   
except Exception as e:
    print(e.args)
    with open(output, "w") as f:
        f.write(str(e.args))
