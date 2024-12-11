from pathlib import Path

import pandas as pd
from snakemake.script import snakemake

# Retrieve parameters and inputs from Snakemake
datasets = snakemake.params.datasets
processed_data_base = snakemake.params.base_path

data = []
for idx, dataset in enumerate(datasets):
    metric_file = Path(
        f"{processed_data_base}/sub-{idx}_{dataset}/ses-0/funcimg/metric.txt"
    )
    error_file = Path(
        f"{processed_data_base}/sub-{idx}_{dataset}/ses-0/funcimg/error.txt"
    )

    # Read metric and error values
    metric = metric_file.read_text().strip() if metric_file.exists() else "N/A"
    error = error_file.read_text().strip() if error_file.exists() else "N/A"

    # Append results
    data.append({"Dataset": dataset, "Metric": metric, "Error": error})

# Create a DataFrame and write to CSV
df = pd.DataFrame(data)
df.to_csv(snakemake.output[0], index=False)
