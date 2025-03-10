from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("calcium-imaging-automation")
except PackageNotFoundError:
    # package is not installed
    pass
