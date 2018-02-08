
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The NCutYX Package

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-NA-6666ff.svg)](https://cran.r-project.org/)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/NCutYX)](https://cran.r-project.org/package=NCutYX)
[![packageversion](https://img.shields.io/badge/Package%20version-0.1.0-orange.svg?style=flat-square)](commits/master)
[![Last-changedate](https://img.shields.io/badge/last%20change-2018--02--08-yellowgreen.svg)](/commits/master)

## Description

The NCutYX package includes functions for clustering genomic data using
graph theory. Each function in this package is a variation on the NCut
measure used to cluster vertices in a graph. The running theme is to use
data sets from different sources and types to improve the clustering
results.

  - The ncut function clusters the columns of a data set using the
    classical normalized cut measure from graph theory.
  - The ancut function clusters one type of data, say gene expressions,
    with the help of a second type of data, like copy number
    aberrations.
  - The muncut function clusters a three-layered graph into K different
    clusters of 3 different data types, say gene expression, copy number
    aberrations and proteins.
  - The pwncut function clusters the columns of X into K clusters by
    giving a weight for each cluster while penalizing them to be similar
    to each other.
  - The mlbncut function works similarly to muncut but it also clusters
    samples into R clusters.
  - The awncut builds similarity matrices for the row of X and an
    assisted dataset Z. Clusters them into K groups while conducting
    feature selection based on the AWNCut method.

To install:

  - latest development version:
    1.  install and load package devtools
    2.  `install_github("Seborinos/NCutYX")`

## NCut

The Normalized Cut (NCut) clusters the columns of Y into K groups using
the NCut graph measure. Builds a similarity matrix for the columns of Y
and clusters them into K groups based on the NCut graph measure.
Correlation, Euclidean and Gaussian distances can be used to construct
the similarity matrix. The NCut measure is minimized using the cross
entropy method, a Monte Carlo optimization technique.

## ANCut

The Assisted NCut (ANcut) clusters the columns of a data set Y into K
groups with the help of an external data set X, which is associated
linearly with Y.

### References:

  - [Hidalgo, Sebastian J. Teran, Mengyun Wu, and Shuangge Ma. “Assisted
    clustering of gene expression data using ANCut.” *BMC genomics* 18.1
    (2017): 623.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5559859/)

## MuNCut

This example shows how to use the muncut function. MuNCut clusters the
columns of data from 3 different sources. It clusters the columns of Z,
Y and X into K clusters by representing each data type as one network
layer. It represents the Z layer depending on Y, and the Y layer
depending on X. Elastic net can be used before the clustering procedure
by using the predictions of Z and Y instead of the actual values to
improve the cluster results. The function muncut will output K clusters
of columns of Z, Y and X.

### References:

  - Sebastian J. Teran Hidalgo and Shuangge Ma. “Clustering Multilayer
    Omics Data using MuNCut.” *Revise and resubmit.*

## PWNCut

The Penalized Weighted NCut (PWNCut) clusters the columns of X into K
clusters by giving a weighted cluster membership while shrinking weights
towards each other.

### References:

  - Sebastian J. Teran Hidalgo, Mengyun Wu and Shuangge Ma. “Penalized
    and weighted clustering of gene expression data using PWNCut.”
    *Submitted.*

## MLBNCut

The Multilayer Biclustering NCut (MLBNCut) clusters the columns and the
rows simultaneously of data from 3 different sources. It clusters the
columns of Z,Y and X into K clusters and the samples into R clusters by
representing each data type as one network layer. It represents the Z
layer depending on Y, and the Y layer depending on X. This function will
output K clusters of columns of Z, Y and X and R clusters of the
samples.

### References:

  - Sebastian J. Teran Hidalgo and Shuangge Ma. “Multilayer Biclustering
    of Omics Data using MLBNCut.” *Work in progress.*

## AWNCut

The Assisted Weighted NCut builds the similarity matrices for the rows
of X and an assisted dataset Z. Clusters them into K groups while
conducting feature selection based on the AWNCut method.

### References:

  - Li, Yang; Bie, Ruofan; Teran Hidalgo, Sebastian; Qin, Yinchen; Wu,
    Mengyun; Ma, Shuangge. “Assisted gene expression-based clustering
    with AWNCut.” *Submitted.*
