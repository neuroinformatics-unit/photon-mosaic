import traceback
from pathlib import Path

from derotation.derotate_batch import derotate
from snakemake.script import snakemake

# Input arguments
read_dataset_path = Path(snakemake.input[0])
output_tif = Path(snakemake.output[0])

output_path_dataset = output_tif.parent
try:
    metrics = derotate(read_dataset_path, output_path_dataset)
    #  save metrics as csv (matrix is already a pandas dataframe)
    metrics.to_csv(output_path_dataset / "derotation_metrics.csv", index=False)

    # make empty error file
    with open(output_path_dataset / "error.txt", "w") as f:
        f.write("")
except Exception:
    with open(output_path_dataset / "error.txt", "w") as f:
        f.write(traceback.format_exc())

    # make empty metrics file
    with open(output_path_dataset / "derotation_metrics.csv", "w") as f:
        f.write("")
