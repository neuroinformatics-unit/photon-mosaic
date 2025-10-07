[![Python Version](https://img.shields.io/pypi/pyversions/photon-mosaic.svg)](https://pypi.org/project/photon-mosaic)
[![PyPI Version](https://img.shields.io/pypi/v/photon-mosaic.svg)](https://pypi.org/project/photon-mosaic)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![CI](https://img.shields.io/github/actions/workflow/status/neuroinformatics-unit/photon-mosaic/test_and_deploy.yml?label=CI)](https://github.com/neuroinformatics-unit/photon-mosaic/actions)
[![Code style: Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/format.json)](https://github.com/astral-sh/ruff)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)


# <img src="docs/source/_static/logo.png" alt="photon-mosaic logo" width="36" style="vertical-align: middle; margin-right:8px"> `photon-mosaic`
`photon-mosaic` is a Snakemake-based toolkit that orchestrates a mosaic of curated algorithms to take you from raw TIFF stacks to downstream, analysis-ready outputs — for example ΔF/F traces, NWB files, or inferred spikes.

It organises algorithms into an automated, self-organising workflow so you can chain preprocessing, registration, signal extraction and post-processing steps into a single, reproducible pipeline. The design focuses on usability for labs that process many imaging sessions and want to scale across an HPC cluster.

<p align="center">
  <img src="https://raw.githubusercontent.com/neuroinformatics-unit/photon-mosaic/refs/heads/improve-docs/docs/source/_static/pm_illustration1.png" alt="photon-mosaic"/>
</p>

## What it does
- Automates multiphoton calcium-imaging analysis: from TIFF stacks to ΔF/F, NWB exports, or spike outputs.
- Arranges a mosaic of community-validated algorithms into reproducible workflows you can configure and extend.
- Runs locally or on an HPC scheduler (SLURM) so you can scale batch processing across sessions.
- Produces standardized outputs and folder layouts that follow the [NeuroBlueprint](https://neuroblueprint.neuroinformatics.dev/latest/index.html) specification.
- Records configs and logs to ensure reproducibility and easy auditing of analysis runs.

## Why use photon-mosaic?
- Self-organising workflows: chain together the algorithms you need and let the pipeline manage execution order and resources.
- Reproducible by design: YAML configs, explicit logs and standardized outputs make it simple to reproduce and share analyses.
- Modular and extensible: integrate alternative preprocessing, deconvolution or denoising modules as needed.
- Built for cluster environments: native support for SLURM lets you process many sessions in parallel.

## Installation

Photon-mosaic requires **Python 3.11** or **3.12**.

```bash
conda create -n photon-mosaic python=3.12
conda activate photon-mosaic
pip install photon-mosaic
```

To install developer tools (e.g., testing and linting):

```bash
pip install 'photon-mosaic[dev]'
```


## Roadmap
### What is implemented today
- Preprocessing: [derotation](https://github.com/neuroinformatics-unit/derotation) and contrast enhancement (see `photon_mosaic/preprocessing`).
- Registration & source extraction: [Suite2p](https://github.com/MouseLand/suite2p).
- Cell detection / anatomical ROI extraction: [Cellpose (v4, including Cellpose-SAM)](https://github.com/MouseLand/cellpose) (Cellpose 4 used by default; Cellpose 3 is also supported).

### Planned additions in the mosaic
- Registration alternative: [NoRMCorre](https://github.com/flatironinstitute/NoRMCorre) implementations for non-rigid motion correction.
- ROI matching: [ROICat](https://github.com/RichieHakim/ROICaT) for inter-session / inter-plane ROI matching.
- Neuropil subtraction / decontamination: methods from the [AllenSDK](https://allensdk.readthedocs.io/en/latest/allensdk.brain_observatory.r_neuropil.html) and [AST-model](https://github.com/znamlab/2p-preprocess).
- Spike deconvolution: [OASIS](https://github.com/j-friedrich/OASIS) and [CASCADE](https://github.com/berenslab/CASCADE) are candidate deconvolution models.

See issues on GitHub for more details and participate in planning.

## Contributing

We welcome issues, feature suggestions, and pull requests. Please refer to our [contribution guidelines](https://photon-mosaic.neuroinformatics.dev/user_guide/index.html) in the documentation for more information.

## References & Links

- [Snakemake Docs](https://snakemake.readthedocs.io/en/stable/)
- [Suite2p Docs](https://suite2p.readthedocs.io/en/latest/)
- [SLURM Executor Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
- [NeuroBlueprint Standard](https://neuroblueprint.neuroinformatics.dev/latest/index.html)
