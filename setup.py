#!/usr/bin/env python

from setuptools import setup

setup(name='s2a',
      version='1.0',
      author='Miles Smith',
      author_email="miles-smith@omrf.org",
      description="Translator to convert Seurat object into AnnData object",
      url='https://gitlab.com/milothepsychic/s2a',
      license='Proprietary',
      classifiers=['Development Status :: 3 - Alpha',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: BSD License',
                   'Programming Language :: Python :: 3.6', 'Programming Language :: Python :: 3.7'],
      packages=["py"],
      keywords=['Seurat', "AnnData", "scRNA-seq"],
      python_requires='>=3.6',
      package_dir={
            's2a': 'py'
      },
      install_requires=['pandas', 'numpy', 'anndata', 'typing'],
      )
