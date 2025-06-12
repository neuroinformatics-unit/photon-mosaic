from importlib.metadata import PackageNotFoundError, version
from importlib.resources import files

try:
    __version__ = version("photon-mosaic")
except PackageNotFoundError:
    # package is not installed
    pass

def get_snakefile_path():
    return files("photon_mosaic").joinpath("workflow", "Snakefile")
