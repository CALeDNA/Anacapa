#!/usr/bin/env Rscript
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("dada2")