# -*- coding: utf-8 -*-
from importlib.metadata import version

"""Top-level package for s2a."""

__author__ = """Miles Smith"""
__email__ = "miles-smith@omrf.org"

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"

from s2a.convert import add_feature_data, add_meta_data, add_reduction
