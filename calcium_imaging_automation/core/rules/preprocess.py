from pathlib import Path

from derotation.analysis.metrics import stability_of_most_detected_blob
from derotation.derotate_batch import derotate
from snakemake.script import snakemake

# Input arguments
read_dataset_path = Path(snakemake.input[0])
output_tif = Path(snakemake.output[0])

output_path_dataset = output_tif.parent.parent
try:
    data = derotate(read_dataset_path, output_path_dataset)
    metric_measured = stability_of_most_detected_blob(data)
    with open(output_path_dataset / "metric.txt", "w") as f:
        f.write(f"dataset: {read_dataset_path.stem} metric: {metric_measured}")
except Exception as e:
    print(e.args)
    with open(output_path_dataset / "error.txt", "w") as f:
        f.write(str(e.args))
