(user_guide/configuration)=
# Configuration

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
