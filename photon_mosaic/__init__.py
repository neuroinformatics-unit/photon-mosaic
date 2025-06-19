from importlib.metadata import PackageNotFoundError, version
from importlib.resources import files

try:
    __version__ = version("photon-mosaic")
except PackageNotFoundError:
    # package is not installed
    pass

def get_snakefile_path():
    """
    Get the path to the Snakemake workflow file.

    Returns
    -------
    Path
        Path to the Snakefile in the photon_mosaic package.
    """
    return files("photon_mosaic").joinpath("workflow", "Snakefile")
