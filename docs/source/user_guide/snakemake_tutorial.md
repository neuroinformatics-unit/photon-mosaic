(user_guide/snakemake_tutorial)=
# Basic Snakemake Tutorial

With `photon-mosaic`, you can run a Snakemake workflow that automatically executes the necessary steps to process your data. The workflow is included in the installed package and can be customized using a YAML configuration file.

The pipeline searches for dataset folders in the specified path and looks for TIFF files in each of them. Each dataset will be processed in parallel, and the results will be saved in a standardized output folder structure under `derivatives`.

## Why use Snakemake?

Snakemake is a powerful workflow management system that allows you to run complex data analysis pipelines in a reproducible and efficient manner. For each defined rule (a rule is a step in the workflow, for instance running Suite2p), Snakemake checks if the output files already exist and whether they are up to date. If not, it runs the rule and generates the outputs.

This approach lets you rerun only the parts of the workflow that need to be updated, avoiding the need to repeat the entire analysis each time.

Here we show examples that do not call directly the `snakemake` command, but instead use the `photon-mosaic` CLI, which is a wrapper around Snakemake that simplifies the execution of the workflow.

## Dry Run
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

## Running the Workflow

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

## Cluster Execution

Once you have tested the workflow locally, you can also submit jobs to a cluster. If you are using SLURM:

```bash
photon-mosaic --jobs 5 --executor slurm
```

## Additional Options

Other useful arguments you can pass:

- `--latency-wait`: wait time before checking if output files are ready
- `--rerun-incomplete`: rerun any incomplete jobs
- `--unlock`: unlock the workflow if it's in a locked state

## Direct Snakemake Usage

You can also run the workflow directly with `snakemake`, using the programmatic path to the bundled Snakefile:

```bash
snakemake --snakefile $(python -c 'import photon_mosaic; print(photon_mosaic.get_snakefile_path())') \
          --configfile path/to/config.yaml \
          --jobs 5
```

This is equivalent to using the `photon-mosaic` CLI but gives full control over the Snakemake interface.

## Path and Wildcard Handling in Snakemake

Key considerations for handling paths and wildcards in Snakemake:

### Wildcard Syntax
Wildcards are used to define the pattern of the output files and are inferred via the outputs paths. Use single curly braces for wildcards in output paths: `{wildcard_name}`. You can use the `wildcards` object to access the values of the wildcards in the rule by constructing a lambda function: `lambda wildcards: str(Path(base_dir) / f"sub-{wildcards.sub_idx}" / "data.npy")`. In the parameters, you can use the wildcards to construct the path to the input files.

Use `pathlib.Path` for cross-platform compatibility and convert paths to strings when using in Snakemake rules.

Example:
```python
rule example:
    input:
        file=lambda wildcards: str(Path(base_dir) / f"sub-{wildcards.sub_idx}" / "data.npy")
    output:
        result=str(Path(output_dir) / "sub-{sub_idx}" / "processed.npy")
    params:
        folder=lambda wildcards: str(Path(base_dir) / f"sub-{wildcards.sub_idx}")
```
