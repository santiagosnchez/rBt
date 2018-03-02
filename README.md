# rBt

## Description

`rBt` includes several functions to handle phylogenetic data. Specifically, it includes a function `read.beast.annot` that can read meta-data from [BEAST1](http://beast.community/index.html) or [BEAST2](https://www.beast2.org) tree files, such as maximum clade credibility (MCC) trees. This package is being developed with phylogenomic species delimitation in mind, but it can be used to many other things like plotting time-trees and visualizing phylogenetic data.

## Installation

Open up an [R](https://www.r-project.org) console:

    install.packages("devtools")
    library(devtools)
    install_github("santiagosnchez/rBt")

## Some functions

* [get.tip.spp](https://github.com/santiagosnchez/rBt/blob/master/R/get.tip.spp.R)
* [mcc.trees2multi](https://github.com/santiagosnchez/rBt/blob/master/R/mcc.trees2multi.R)
* [mcc2](https://github.com/santiagosnchez/rBt/blob/master/R/mcc2.R)
* [order.edges](https://github.com/santiagosnchez/rBt/blob/master/R/order.edges.R)
* [read.beast.annot](https://github.com/santiagosnchez/rBt/blob/master/R/read.beast.annot.R)
* [which.node](https://github.com/santiagosnchez/rBt/blob/master/R/which.node.R)
* [plot.phylo.HPD](https://github.com/santiagosnchez/rBt/blob/master/R/plot.phylo.HPD.R)

## Examples

Run the function with `?` in front (e.g. `?mcc2`) to access the help files, which have examples and more information.

## Citation

Sánchez-Ramírez, S. (2018) rBt: Handy functions for dealing with BEAST trees in R, GitHub repository: https://github.com/santiagosnchez/rBt
