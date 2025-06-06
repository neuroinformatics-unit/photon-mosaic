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

Photon-mosaic requires **Python 3.11** or **3.12** and depends on a custom fork of Suite2p for compatibility.

```bash
conda create -n photon-mosaic python=3.12
conda activate photon-mosaic
pip install photon-mosaic
pip install git+https://github.com/neuroinformatics-unit/suite2p.git
```
N.B.: as you can see, we are using a custom fork of Suite2p to ensure compatibility with the latest Python versions and to include additional features. This fork is maintained by the Neuroinformatics Unit.

To install developer tools (e.g., testing and linting):

```bash
pip install -r requirements-dev.txt
```

## Configuration

On first run, photon-mosaic will create a user config at `~/.photon_mosaic/config.yaml` if it does not exist.

You can:
- Edit this file directly for user-wide defaults.
- Override key paths at the command line:
  ```bash
  photon-mosaic --raw_data_base /my/data --processed_data_base /my/processed --jobs 5
  ```
- Use a project-specific config:
  ```bash
  photon-mosaic --config ./my/path/to/config.yaml --jobs 5
  ```

**Note:**
- The config used for the run (with any overrides) is always exported to a timestamped file in the `derivatives/photon-mosaic/configs/` directory.
- Snakemake logs are always dumped to a timestamped file in the `derivatives/photon-mosaic/logs/` directory.
- Both logs and configs are organized with timestamps (format: YYYYMMDD_HHMMSS) for easy tracking of different runs.

Here is an example of a `config.yaml` file:

```yaml
raw_data_base: "/path/to/raw/"
processed_data_base: "/path/to/processed/"

suite2p_ops:
  fs: 7.5
  nplanes: 1
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

If you don't have access to a cluster or SLURM, set `use_slurm: false` to run locally.

## Basic snakemake tutorial

With `photon-mosaic`, you can run a Snakemake workflow that automatically executes the necessary steps to process your data. The workflow is included in the installed package and can be customized using a YAML configuration file.

The pipeline searches for dataset folders in the specified path and looks for TIFF files in each of them. Each dataset will be processed in parallel, and the results will be saved in a standardized output folder structure under `derivatives`.

### Why use Snakemake?

Snakemake is a powerful workflow management system that allows you to run complex data analysis pipelines in a reproducible and efficient manner. For each defined rule (a rule is a step in the workflow, for instance running Suite2p), Snakemake checks if the output files already exist and whether they are up to date. If not, it runs the rule and generates the outputs.

This approach lets you rerun only the parts of the workflow that need to be updated, avoiding the need to repeat the entire analysis each time.

Here we show examples that do not call directly the `snakemake` command, but instead use the `photon-mosaic` CLI, which is a wrapper around Snakemake that simplifies the execution of the workflow.

**Dry Run**
A dry run is a simulation that shows what would happen if the workflow were executed, without actually running any commands. This is useful for verifying that everything is set up correctly. The output includes a DAG (directed acyclic graph) showing dependencies between rules, which files will be created, and which rules will be executed.

To preview the workflow without running it:

```bash
photon-mosaic --jobs 1 --dry-run
```

By default, the workflow uses the configuration file included in the package. To run with your own configuration:

```bash
photon-mosaic --config path/to/config.yaml --jobs 1 --dry-run
```

`--jobs` specifies the number of jobs to run in parallel. You can increase this number to parallelize execution across datasets. A dry run can also be abbreviated to `-np` if using Snakemake directly.

To run the full workflow:

```bash
photon-mosaic --jobs 5
```

To force the re-execution of a specific rule:

```bash
photon-mosaic --jobs 5 --forcerun suite2p
```

To reprocess a specific dataset, you can specify a target output file (e.g., `F.npy`):

```bash
photon-mosaic --jobs 1 /path/to/derivatives/dataset_name/suite2p/plane_0/F.npy
```

Once you have tested the workflow locally, you can also submit jobs to a cluster. If you are using SLURM:

```bash
photon-mosaic --jobs 5 --executor slurm
```

Other useful arguments you can pass:

- `--latency-wait`: wait time before checking if output files are ready
- `--rerun-incomplete`: rerun any incomplete jobs
- `--unlock`: unlock the workflow if it's in a locked state

You can also run the workflow directly with `snakemake`, using the programmatic path to the bundled Snakefile:

```bash
snakemake --snakefile $(python -c 'import photon_mosaic; print(photon_mosaic.get_snakefile_path())') \
          --configfile path/to/config.yaml \
          --jobs 5
```

This is equivalent to using the `photon-mosaic` CLI but gives full control over the Snakemake interface.


## Contributing

We welcome issues, feature suggestions, and pull requests. Please refer to our [contribution guidelines](https://photon-mosaic.neuroinformatics.dev/user_guide/index.html) in the documentation for more information.

## References & Links

- [Snakemake Docs](https://snakemake.readthedocs.io/en/stable/)
- [Suite2p Docs](https://suite2p.readthedocs.io/en/latest/)
- [Custom Suite2p Fork](https://github.com/neuroinformatics-unit/suite2p.git)
- [SLURM Executor Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
- [NeuroBlueprint Standard](https://neuroblueprint.neuroinformatics.dev/latest/index.html)
