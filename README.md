# MetroSCREEN

![PyPI](https://img.shields.io/pypi/v/MetroSCREEN)
![Downloads](https://pepy.tech/badge/MetroSCREEN)
![Documentation Status](https://readthedocs.org/projects/MetroSCREEN/badge/?version=latest)

Cell identification in high-resolution Spatial Transcriptomics

MetroSCREEN is a computational method to perform cell segmentation on high-resolution spatial transcriptomics (ST) data, including sequncing-based (e.g. Stereo-seq and Seq-Scope) and imaging-based (e.g. seqFISH+ and STARmap) technologies.

![avatar](docs/_static/img/MetroSCREEN_framework.png)

## Change Log
### v0.0.1a
* Build MetroSCREEN.
### v1.0.0
* Release MetroSCREEN.


## Install MetroSCREEN
```bash
devtools::install_github('wanglabtongji/MetroSCREEN')
```
## Installation of Other Dependencies
Install the osqp package for optimization using install.packages('osqp'), If you encounter any issue during MetroSCREEN installation.
Install the dplyr package using install.packages('dplyr').
For single cell data analysis, we provide pipeline to work with Seurat(Seurat V4). Please install Seurat package by install.packages('Seurat').

## Documentation
For full installation and usage of MetroSCREEN, please refer to the [documentation](https://metroscreen-rtd-rutorial.readthedocs.io/en/latest/).
