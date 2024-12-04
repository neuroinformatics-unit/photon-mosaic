from pathlib import Path

from derotation.analysis.metrics import stability_of_most_detected_blob
from derotation.derotate_batch import derotate
from snakemake.script import snakemake

try:
    # Input arguments
    read_dataset_path = Path(snakemake.input[0])
    write_dataset_path = Path(snakemake.input[1])
    output = snakemake.output[0]

    output_path_dataset = write_dataset_path / "ses-0/funcimg/"

    data = derotate(read_dataset_path, output_path_dataset)
    metric_measured = stability_of_most_detected_blob(data)
    with open(output, "w") as f:
        f.write(f"dataset: {read_dataset_path.stem} metric: {metric_measured}")
except Exception as e:
    print(e.args)
    with open(output, "w") as f:
        f.write(str(e.args))
