# photon-mosaic

**photon-mosaic** is a package in development for the automated analysis of multi-photon calcium imaging datasets. It provides a toolbox of commonly used processing and analysis methods, allowing users to create their own processing pipelines.

## Installation

To install **photon-mosaic**, clone the repository, create a virtual environment (e.g., using Conda), and install the package:

```bash
git clone https://github.com/neuroinformatics-unit/photon-mosaic.git
cd photon-mosaic
conda create -n photon-mosaic python=3.12
conda activate photon-mosaic
pip install .
```

### Suite2p Dependency

Currently, **photon-mosaic** requires a custom fork of the `suite2p` package. Install it using:

```bash
pip install git+https://github.com/neuroinformatics-unit/suite2p.git
```

This is necessary because the official `suite2p` package does not yet support the latest versions of Python and NumPy.

## Usage

This package leverages **Snakemake** to create analysis pipelines. A minimal example pipeline is provided to process multiple datasets in parallel using `suite2p`, taking advantage of local cluster resources via the Slurm plugin for Snakemake.

### Configuring the Pipeline

To run the pipeline, edit the **Snakemake file** (`Snakefile`) located in the `workflow/` folder. You need to specify:

- The **raw data folder**, which should contain subfolders named after subjects.
  - Each subject folder should contain one `.tif` file per experimental session.
- The **output folder** where processed data will be saved.
- The **Suite2p options file** location.
- The **SLURM resources**, including partition, number of cores, and memory allocation.

The output data will be formatted following the **NeuroBlueprint standard**.

### Running the Pipeline

To execute the pipeline, run:

```bash
snakemake --executor slurm --jobs 20 --latency-wait 10 all
```

- `--jobs 20`: Specifies the number of subjects to process in parallel.
- `--latency-wait 10`: Waits 10 seconds for output files to be created.
- `all`: Runs the full pipeline defined in the `Snakefile`.

#### Dry Run Mode

To preview the execution plan without running any processing, use:

```bash
snakemake --executor slurm --jobs 20 --latency-wait 10 --dry-run all
```

#### Customizing Suite2p Options

Modify the Suite2p configuration in:

```
photon_mosaic/rules/s2p_options.py
```

## Contributing

Contributions are welcome! If you have suggestions or questions, please open an issue on the repository.

## References

- [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/)
- [Suite2p Documentation](https://suite2p.readthedocs.io/en/latest/)
- [Custom Suite2p Fork](https://github.com/neuroinformatics-unit/suite2p.git)
- [Slurm for Snakemake](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
- [NeuroBlueprint Standard](https://neuroblueprint.neuroinformatics.dev/latest/index.html)

