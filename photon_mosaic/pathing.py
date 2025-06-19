import os


def cross_platform_path(path):
    """
    Convert path to string format appropriate for the current platform.
    On Windows, uses forward slashes (as_posix()) for Snakemake compatibility.
    On Unix-like systems, uses native path separators (str()).
    """
    if os.name == "nt":
        return path.as_posix()
    else:
        return str(path)
