# rBt

## Description

`rBt` includes several functions to handle phylogenetic data. Specifically, it includes a function `read.beast.annot` that can read meta-data from [BEAST1](http://beast.community/index.html) or [BEAST2](https://www.beast2.org) tree files, such as maximum clade credibility (MCC) trees. This package is being developed with phylogenomic species delimitation in mind.

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


