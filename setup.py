#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

with open("requirements.txt") as requirements_file:
    requirements = requirements_file.read().splitlines()

setup_requirements = []

test_requirements = []

setup(
    author="Miles Smith",
    author_email="miles-smith@omrf.org",
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    description="glue to link Seurat to anndata objects",
    install_requires=requirements,
    license="BSD license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords=["Seurat", "AnnData", "scRNA-seq", "Scanpy"],
    name="s2a",
    packages=find_packages(include=["s2a", "s2a.*"]),
    package_dir={"s2a": "s2a"},
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://gitlab.com/milothepsychic/s2a",
    version="0.1.0",
    zip_safe=False,
)
