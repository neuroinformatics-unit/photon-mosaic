# `photon-mosaic`

`photon-mosaic` is a Snakemake-based pipeline for the automated and reproducible analysis of multiphoton calcium imaging datasets. It currently integrates [Suite2p](https://suite2p.readthedocs.io/en/latest/) for image registration and signal extraction, with plans to support additional analysis modules in the future.

<p align="center">
  <img src="https://raw.githubusercontent.com/neuroinformatics-unit/photon-mosaic/main/docs/source/_static/photon-mosaic.png" alt="photon-mosaic" width="30%"/>
</p>

&nbsp;
## Overview
`photon-mosaic` can leverage SLURM job scheduling, allows standardized and reproducible workflows configurable via a simple YAML config file and produces standardized output folder structures following the [NeuroBlueprint](https://neuroblueprint.neuroinformatics.dev/latest/index.html) specification.

This tool is especially suited for labs that store their data on servers directly connected to an HPC cluster and want to batch-process multiple imaging sessions in parallel.

The current structure sets the stage for future modular integration of preprocessing, neuropil decontamination and deconvolution of choice, and more.

## Installation

Photon-mosaic requires **Python 3.11** or **3.12** and installs a custom fork of Suite2p for compatibility.

```bash
git clone https://github.com/neuroinformatics-unit/photon-mosaic.git
cd photon-mosaic
conda create -n photon-mosaic python=3.12
conda activate photon-mosaic
pip install .
```

To install developer tools (e.g., testing and linting):

```bash
pip install '.[dev]'

## Configuration

Edit or create your `config.yaml` file like so:

```yaml
raw_data_base: "/path/to/raw/"
processed_data_base: "/path/to/processed/"

suite2p_ops:
  fs: 7.5
  nplanes: 2
  tau: 0.8
  nonrigid: true
  diameter: 10

slurm:
  use_slurm: true
  partition: "cpu"
  mem_mb: 16000
  cpu_per_task: 1
  tasks: 1
  nodes: 1
```

If you don’t have access to a cluster or SLURM, set `use_slurm: false` to run locally.

## Basic snakemake tutorial

With `snakemake` you can run a workflow that automatically runs the necessary steps to process your data. The workflow is pre-defined in `workflow/Snakefile` and can be customized using the provided configuration file.

It searches for dataset folders in the specified path and searches for tiffs in each of them. Each dataset will be processed in parallel and the results will be saved in the specified output folder called `derivatives`.

### Why use snakemake?
Snakemake is a powerful workflow management system that allows you to run complex data analysis pipelines in a reproducible and efficient manner. For each defined rule (a rule is a step in the workflow, for instance running suite2p), Snakemake will check if the output files already exist and if they are up to date. If they are not, it will run the rule and create the output files. This way, you can easily rerun only the parts of the workflow that need to be updated, without having to rerun the entire analysis pipeline each time.

**Dry Run**
A dry run is a simulation of the workflow that shows you what would happen if you ran it, without actually executing any commands. This is useful for checking if everything is set up correctly before running the workflow. What you will see as an output is a DAG, i.e. a directed acyclic graph, that shows the dependencies between the different rules in the workflow. You can also see which files will be created and which rules will be executed. For rules to be linked together, input and output names must match: rule A will be linked to rule B if the output of rule A is the input of rule B.

To preview the workflow without running it:
```bash
snakemake --jobs 1 all --dry-run
```
`all` is a keyword that tells snakemake to run all the rules in the workflow.
The `--jobs` argument specifies the number of jobs to run in parallel. In this case, we are running one job at a time. You can increase this number to run multiple jobs in parallel.
Dry run can also be abbreviated to `-np`.

To run the workflow you can skip the `--dry-run` argument and run the command:
```bash
snakemake --jobs 5 all
```

If you wish to force the re-execution of a given rule you can use the `--force-rerun` argument followed by the name of the rule you want to rerun. For example, if you want to rerun the rule `suite2p`, you can use the command:
```bash
snakemake --jobs 5 all --forcerun suite2p
```

You can also rerun a specific dataset by specifying the output of interest:
```bash
snakemake --jobs 1 /path/to/derivatives/dataset_name/suite2p/plane_0/F.npy
```
This will trigger the analysis that lead to the creation of the file `F.npy` in the specified dataset folder.

Once you have tested the workflow locally, you can use further arguments to submit the jobs to a cluster. If you are using SLURM, you can use the following command:

```bash
snakemake --executor slurm --jobs 5 all
```

Other useful arguments are:
- `--latency-wait`: to wait for a certain amount of time before checking if the output files are ready.
- `--rerun-incomplete`: to rerun any incomplete jobs.
- `--unlock`: to unlock the workflow.

## Contributing

We welcome issues, feature suggestions, and pull requests. Please refer to our [contribution guidelines](https://photon-mosaic.neuroinformatics.dev/user_guide/index.html) in the documentation for more information.

## References & Links

- [Snakemake Docs](https://snakemake.readthedocs.io/en/stable/)
- [Suite2p Docs](https://suite2p.readthedocs.io/en/latest/)
- [Custom Suite2p Fork](https://github.com/neuroinformatics-unit/suite2p.git)
- [SLURM Executor Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
- [NeuroBlueprint Standard](https://neuroblueprint.neuroinformatics.dev/latest/index.html)
