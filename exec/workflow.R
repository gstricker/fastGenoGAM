## library(GenomicRanges)
## library(SummarizedExperiment)
## library(Matrix)
## library(data.table)
## library(futile.logger)
## library(profvis)

## source("./R/GenoGAMSettings-class.R")
## source("./R/GenoGAMDataSet-class.R")
## source("./R/GenoGAMSetup-class.R")
## source("./R/readData.R")
## source("./R/sf.R")
## source("./R/cv.R")

library(fastGenoGAM)
folder <- "/s/project/coreProm/Michi/thornton/align_STAR"
config <- "/data/ouga/home/ag_gagneur/strickeg/workspace/analysis/diffBinding/config.txt"

## folder <- "/s/project/coreProm/data/sacCer2"
## config <- "/data/ouga/home/ag_gagneur/strickeg/workspace/analysis/peakcalls/config_tfiib.txt"

settings <- GenoGAMSettings()
slot(settings, "optimControl")$betaMaxit <- 10000
ggd <- GenoGAMDataSet(config, 10000, 250, design = ~ s(x) + s(x, by = genotype), directory = folder, settings = settings)
ggd <- computeSizeFactors(ggd)

lambda <- 9030.333292
theta <- 4.670743
##debug(genogam)
result <- genogam(ggd, intervalSize = 200, lambda = lambda, theta = theta)
