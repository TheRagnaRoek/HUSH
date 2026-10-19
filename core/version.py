"""The installed HUSH version, read from the package metadata (pyproject.toml)."""

from importlib import metadata

try:
    __version__ = metadata.version("hush")
except metadata.PackageNotFoundError:  # pragma: no cover - not installed
    __version__ = "0.0.0"
