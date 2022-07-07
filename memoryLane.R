## Bukola Ajanaku
## July 26, 2021
## memoryLane.R: for alzheimer's DNA methylation data
## on screen 2 as toRemember
## module load udunits proj gdal geos
## module load R/4.0.2
## R
## save(list = ls(), file = "/sc/arion/work/ajanab01/newFinalPlayData.RDATA")

.libPaths(c("/hpc/users/ajanab01/.Rlib", .libPaths()))

suppressPackageStartupMessages(library(decorate))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(BiocParallel))

options(xtable.type="html")

knitr::opts_chunk$set(
  echo=TRUE,
  warning=FALSE,
  message=TRUE,
  error = FALSE,
  tidy = FALSE,
  cache = TRUE,
  cache.lazy = FALSE,
  dev = c("png", "pdf"),
  fig.width=7, fig.height=7)

options(markdown.HTML.stylesheet = 'css/custom.css')

featureLocation <- read.csv("/sc/arion/projects/epigenAD/Bukola/H3K9ac/H3K9acDomains.csv", sep="\t", header=TRUE)
featureLocation <- as(featureLocation, "GRanges")

response = read.csv("/sc/arion/projects/epigenAD/Bukola/H3K9ac/ReadCounts.csv", sep="\t", header=TRUE)

metadata = read.table("/sc/arion/projects/epigenAD/Bukola/H3K9ac/SYNAPSE_METADATA_MANIFEST.tsv", sep = '\t', header = TRUE)
