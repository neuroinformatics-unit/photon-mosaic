[![Python Version](https://img.shields.io/pypi/pyversions/photon-mosaic.svg)](https://pypi.org/project/photon-mosaic)
[![PyPI Version](https://img.shields.io/pypi/v/photon-mosaic.svg)](https://pypi.org/project/photon-mosaic)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![CI](https://img.shields.io/github/actions/workflow/status/neuroinformatics-unit/photon-mosaic/test_and_deploy.yml?label=CI)](https://github.com/neuroinformatics-unit/photon-mosaic/actions)
[![Code style: Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/format.json)](https://github.com/astral-sh/ruff)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)


# <img src="docs/source/_static/logo.png" alt="photon-mosaic logo" width="36" style="vertical-align: middle; margin-right:8px"> `photon-mosaic`
`photon-mosaic` is a Snakemake-based toolkit for processing multiphoton datasets. It orchestrates a curated collection of algorithms to transform your raw data (e.g., TIFF files) into analysis-ready outputs, such as ΔF/F traces, NWB files, or inferred spikes.

Each analysis step is integrated into an automated workflow, allowing you to chain preprocessing, registration, signal extraction, and post-processing steps into a single, reproducible pipeline. The design prioritizes usability for labs that process many imaging sessions and need to scale across an HPC cluster.

<p align="center">
  <img src="https://raw.githubusercontent.com/neuroinformatics-unit/photon-mosaic/refs/heads/improve-docs/docs/source/_static/pm_illustration1.png" alt="photon-mosaic"/>
</p>

This is made possible by [Snakemake](https://snakemake.readthedocs.io/en/stable/), a workflow management system that provides a powerful and flexible framework for defining and executing complex data processing pipelines. Snakemake automatically builds a directed acyclic graph (DAG) of all the steps in your analysis, ensuring that each step is executed in the correct order and that intermediate results are cached to avoid redundant computations. `photon-mosaic` also includes a [SLURM executor plugin for Snakemake](https://github.com/snakemake/snakemake-executor-plugin-slurm) to seamlessly scale your analysis across an HPC cluster. To ensure consistency and reproducibility, `photon-mosaic` writes processed data according to the [NeuroBlueprint](https://neuroblueprint.neuroinformatics.dev/latest/index.html) data standard for organizing and storing multiphoton imaging data. 

The goal of `photon-mosaic` is to provide a modular, extensible, and user-friendly framework for multiphoton data analysis that can be easily adapted to different experimental designs and analysis requirements. For each processing step, we aim to vet and integrate the best available open-source tools, providing sensible defaults tailored to the specific data type and experimental modality.

## Roadmap
### Current features
- **Preprocessing**: [derotation](https://github.com/neuroinformatics-unit/derotation) and contrast enhancement (see `photon_mosaic/preprocessing`).
- **Registration & source extraction**: [Suite2p](https://github.com/MouseLand/suite2p).
- **Cell detection / anatomical ROI extraction**: [Cellpose (v3 or v4, including Cellpose-SAM)](https://github.com/MouseLand/cellpose).

### Planned additions
- **Registration alternative**: [NoRMCorre](https://github.com/flatironinstitute/NoRMCorre) implementations for non-rigid motion correction.
- **ROI matching**: [ROICat](https://github.com/RichieHakim/ROICaT) for inter-session / inter-plane ROI matching.
- **Neuropil subtraction / decontamination**: Methods from the [AllenSDK](https://allensdk.readthedocs.io/en/latest/allensdk.brain_observatory.r_neuropil.html) and [AST-model](https://github.com/znamlab/2p-preprocess).
- **Spike deconvolution**: [OASIS](https://github.com/j-friedrich/OASIS) and [CASCADE](https://github.com/HelmchenLabSoftware/Cascade).

## Installation

Photon-mosaic requires **Python 3.11** or **3.12**.

```bash
conda create -n photon-mosaic python=3.12
conda activate photon-mosaic
pip install photon-mosaic
```

To install with developer tools (e.g., testing and linting):

```bash
pip install 'photon-mosaic[dev]'
```

## Contributing

We welcome issues, feature suggestions, and pull requests. Please refer to our [contribution guidelines](https://photon-mosaic.neuroinformatics.dev/contributing.html) in the documentation for more details.
