# rBt

## Description

`rBt` includes several functions to handle phylogenetic data. Specifically, it includes a function `read.beast.annot` that can read meta-data from [BEAST1](http://beast.community/index.html) or [BEAST2](https://www.beast2.org) tree files, such as maximum clade credibility (MCC) trees. This package is being developed with phylogenomic species delimitation in mind, but it can be used to many other things like plotting time-trees and visualizing phylogenetic data.

## Installation

Open up an [R](https://www.r-project.org) console:

    install.packages("devtools")
    library(devtools)
    install_github("santiagosnchez/rBt")

## Functions

* get.tip.spp
* mcc.trees2multi
* mcc2
* order.edges
* read.beast.annot
* which.node
* plot.phylo.HPD

## Examples

Run the function with `?` in front (e.g. `?mcc2`) to access the help files, which have examples and more information.

## Citation

Sánchez-Ramírez, S. (2018) rBt: Handy functions for dealing with BEAST trees in R, GitHub repository: https://github.com/santiagosnchez/rBt
