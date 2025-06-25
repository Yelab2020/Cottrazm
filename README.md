# Cottrazm: Construct Tumor Transition Zone (Boundary) Microenvironment based on spatial transcriptomics
## Overview

Cottrazm (**Co**nstruct **T**umor **Tra**nsition **Z**one (Boundary) **M**icroenvironment based on spatial transcriptomics) aims to construct the microenvironment of tumor boundary based on spatial transcriptomics, single-cell transcriptomics and HE-stained histological images. It consists of three core functions: determining the tumor boundary (Cottrazm-BoundaryDefine), deconvoluting spatial transcriptomics (Cottrazm-SpatialDecon), and reconstructing a spatial gene expression matrix for sub-spots (Cottrazm-SpatialRecon).

 Taken together, Cottrazm provides an integrated tool framework to dissect the tumor spatial microenvironment and facilitates the discovery of functional biological insights, thereby identifying therapeutic targets in oncologic ST datasets.
 

 ![overview_image](https://github.com/Yelab2020/Cottrazm/blob/main/doc/github%20figure1.png)
## Installation

The development version can be installed from git hub with 'devtools'.
```
# Install devtools, if necessary
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("Yelab2020/Cottrazm")
```
You can also download Cottrazm through docker. Make sure you have docker installed.

```
docker pull skye12138/cottrazm:latest
```

## Usage

For usage of Cottrazm, please refer to the [vignette](doc/my-vignette.pdf)

## How to cite Cottrazm

Please cite the following manuscript:
Xun, Z., Ding, X., Zhang, Y. et al. Reconstruction of the tumor spatial microenvironment along the malignant-boundary-nonmalignant axis. Nat Commun 14, 933 (2023). https://doi.org/10.1038/s41467-023-36560-7
