# calcium-imaging-automation

CIMAT simplifies the analysis of multi-photon calcium imaging data by integrating algorithms from tools like Suite2p and Caiman into a modular Snakemake pipeline. Researchers can evaluate, compare, and combine methods for each processing step, such as registration or source extraction, and explore metrics to identify the best fit for their datasets.

With support for local or cluster-based parallelization, CIMAT provides visualization tools, reports, and guides to streamline decision-making and enhance reproducibility.

## Installation


### Run workflow with Snakemake
To extract dataset names
```bash
snakemake --cores 1 setup_output.txt
```

Run all jobs in the pipeline:
```bash
snakemake --executor slurm --jobs 20 --latency-wait 10 all
```
Add `-np --printshellcmds` for a dry run with commands printed to the terminal.
