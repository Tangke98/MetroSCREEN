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
git clone https://github.com/wanglabtongji/MetroSCREEN.git
cd MetroSCREEN
conda create -n MetroSCREEN python=3.10
pip install -r requirements.txt
pip install .
```

## Documentation
For full installation and usage of MetroSCREEN, please refer to the [documentation](https://MetroSCREEN.readthedocs.io/en/latest/).

## Usage
```bash
MetroSCREEN --help
usage: MetroSCREEN [-h] [-v] {seg,align,watershed,impute} ...

MetroSCREEN (Cell identification in high-resolution Spatial Transcriptomics) is a cell segmentation tool for high-resolution spatial transcriptomics.

positional arguments:
  {seg,align,watershed,impute}
    seg                 Run MetroSCREEN segmentation on high-resolution spatial transcriptomics.
    align               Refine alignment between image and spatial transcriptomics.
    watershed           Run initial watershed segmentation on the staining image.
    impute              Perform spatially-aware gene imputation within each cluster.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version info.
```
