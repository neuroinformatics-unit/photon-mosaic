# calcium-imaging-automation

CIMAT simplifies the analysis of multi-photon calcium imaging data by integrating algorithms from tools like Suite2p and Caiman into a modular Snakemake pipeline. Researchers can evaluate, compare, and combine methods for each processing step, such as registration or source extraction, and explore metrics to identify the best fit for their datasets.

With support for local or cluster-based parallelization, CIMAT provides visualization tools, reports, and guides to streamline decision-making and enhance reproducibility.

## Installation


### Run workflow with Snakemake
Run all jobs in the pipeline:
```bash
snakemake --executor slurm --jobs 20 --latency-wait 10 all --forcerun preprocess --rerun-incomplete
```
Add `-np --printshellcmds` for a dry run with commands printed to the terminal.

### See interactive report with datavzrd
Build the csv:
```bash
snakemake --cores 1 workflow/results/data/summary.csv
```
Create the report:
```bash
datavzrd workflow/resources/datavzrd_config.yaml --output workflow/results/datavzrd
```
Then open the report (`index.html`) in a browser.

specific dataset, rerun:
```bash
snakemake --executor slurm --jobs 20 --latency-wait 10     /ceph/margrie/laura/cimaut/derivatives/sub-1_230802CAA1120182/ses-0/funcimg/derotation/derotated_full.tif     --forcerun preprocess --rerun-incomplete
```

summary plot:
```bash
snakemake --cores 1 --latency-wait 10 workflow/results/data/stability_metric.png
```
