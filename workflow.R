library(GenomicRanges)
library(SummarizedExperiment)
library(Matrix)
library(data.table)
library(futile.logger)
library(profvis)

source("./R/GenoGAMSettings-class.R")
source("./R/GenoGAMDataSet-class.R")
source("./R/GenoGAMSetup-class.R")
source("./R/readData.R")
source("./R/sf.R")
source("./R/cv.R")

flog.threshold(DEBUG)
folder <- "/s/project/coreProm/Michi/thornton/align_STAR"
config <- "/data/ouga/home/ag_gagneur/strickeg/workspace/analysis/diffBinding/config.txt"

ggd <- GenoGAMDataSet(config, 2000, 250, design = ~ s(x), directory = folder)
ggd <- computeSizeFactors(ggd)

