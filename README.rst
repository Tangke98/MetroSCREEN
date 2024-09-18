# Cellist

![PyPI](https://img.shields.io/pypi/v/cellist)
![Downloads](https://pepy.tech/badge/cellist)
![Documentation Status](https://readthedocs.org/projects/cellist/badge/?version=latest)

Cell identification in high-resolution Spatial Transcriptomics

Cellist is a computational method to perform cell segmentation on high-resolution spatial transcriptomics (ST) data, including sequncing-based (e.g. Stereo-seq and Seq-Scope) and imaging-based (e.g. seqFISH+ and STARmap) technologies.

![avatar](docs/_static/img/Cellist_workflow.png)

## Change Log
### v0.0.1a
* Build Cellist.
### v1.0.0
* Release Cellist.


## Install Cellist
```bash
git clone https://github.com/wanglabtongji/Cellist.git
cd Cellist
conda create -n Cellist python=3.10
pip install -r requirements.txt
pip install .
```

## Documentation
For full installation and usage of Cellist, please refer to the [documentation](https://cellist.readthedocs.io/en/latest/).

## Usage
```bash
cellist --help
usage: cellist [-h] [-v] {seg,align,watershed,impute} ...

Cellist (Cell identification in high-resolution Spatial Transcriptomics) is a cell segmentation tool for high-resolution spatial transcriptomics.

positional arguments:
  {seg,align,watershed,impute}
    seg                 Run Cellist segmentation on high-resolution spatial transcriptomics.
    align               Refine alignment between image and spatial transcriptomics.
    watershed           Run initial watershed segmentation on the staining image.
    impute              Perform spatially-aware gene imputation within each cluster.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version info.
```
