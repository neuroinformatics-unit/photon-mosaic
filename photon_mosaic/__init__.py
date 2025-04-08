from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("photon-mosaic")
except PackageNotFoundError:
    # package is not installed
    pass
