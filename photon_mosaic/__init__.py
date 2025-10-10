from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("photon-mosaic")
except PackageNotFoundError:
    # package is not installed
    pass

# Import snakemake utilities for convenience
from photon_mosaic.snakemake_utils import (
    cross_platform_path,
    get_snakefile_path,
    log_cuda_availability,
)

__all__ = [
    "get_snakefile_path",
    "cross_platform_path",
    "log_cuda_availability",
]
