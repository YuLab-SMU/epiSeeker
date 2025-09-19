# epiSeeker: an R package for Annotation, Comparison and Visualization of multi-omics epigenetic data

<!-- r badge_bioc_release("epiSeeker", "green")` 
r badge_devel("YuLab-SMU/epiSeeker", "green")`
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/epiSeeker.svg)](https://www.bioconductor.org/packages/devel/bioc/html/epiSeeker.html#since)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
&#10;
[![platform](http://www.bioconductor.org/shields/availability/devel/epiSeeker.svg)](https://www.bioconductor.org/packages/devel/bioc/html/epiSeeker.html#archives)
[![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/epiSeeker.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/epiSeeker/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2025--09--19-green.svg)](https://github.com/YuLab-SMU/epiSeeker/commits/master)
-->

This package implements functions to analyzing multi-omics epigenetic
data. Data of fragment type and base type are supported by epiSeeker. It
provide function to retrieve the nearest genes around the peak, annotate
genomic region of the peak, statstical methods for estimate the
significance of overlap among peak data sets, and motif analysis. It
incorporate GEO database for user to compare the own dataset with those
deposited in database. The comparison can be used to infer cooperative
regulation and thus can be used to generate hypotheses. Several
visualization functions are implemented to summarize the coverage of the
peak experiment, average profile and heatmap of peaks binding to TSS
regions, genomic annotation, distance to TSS, overlap of peaks or genes,
and the single-base resolution epigenetic data by considering the
strand, motif, and additional information.

## :writing_hand: Authors

Guangchuang YU

School of Basic Medical Sciences, Southern Medical University

<https://yulab-smu.top>

## :arrow_double_down: Installation

Get the released version from Bioconductor:

``` r
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("epiSeeker")
```

Or the development version from github:

``` r
## install.packages("devtools")
devtools::install_github("YuLab-SMU/epiSeeker")
```

## Contributing

We welcome any contributions! By participating in this project you agree
to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
