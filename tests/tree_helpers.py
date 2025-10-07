"""Small helpers for printing directory trees in tests.

Provides a compact generator-style `tree()` implementation that yields a
visual tree (lines prefixed with ├──, └── and indented with │) for a
given directory. This is intentionally lightweight and has no external
dependencies so it can be imported from test code without extra setup.
"""

from pathlib import Path
from typing import Generator

# prefix components:
space = "    "
branch = "│   "
# pointers:
tee = "├── "
last = "└── "


def tree(dir_path: Path, prefix: str = "") -> Generator[str, None, None]:
    """Yield a visual tree structure for ``dir_path`` line by line.

    Each yielded line is prefixed by ``prefix`` so callers can accumulate
    indent levels when recursing. Entries are sorted by name for
    determinism in tests.
    """
    try:
        contents = sorted(list(dir_path.iterdir()), key=lambda p: p.name)
    except Exception:
        return

    # contents each get pointers that are ├── with a final └── :
    if not contents:
        return

    pointers = [tee] * (len(contents) - 1) + [last]
    for pointer, path in zip(pointers, contents):
        yield prefix + pointer + path.name
        if path.is_dir():  # extend the prefix and recurse:
            extension = branch if pointer == tee else space
            # i.e. space because last, └── , above so no more |
            yield from tree(path, prefix=prefix + extension)
