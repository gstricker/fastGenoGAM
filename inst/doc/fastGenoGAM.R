## ----setup, echo=FALSE, results="hide"-------------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png", 
                      message=FALSE, error=FALSE, warning=FALSE)

## --------------------------------------------------------------------------
library(fastGenoGAM)

## A.
## specify folder and experiment design path
wd <- system.file("extdata/Set1", package='fastGenoGAM')
folder <- file.path(wd, "bam")
expDesign <- file.path(wd, "experimentDesign.txt")

## B. 
## set HDF5 folder where the data will be stored
## Note, don't usually use /tmp because all your data and
## results get deleted otherwise later
hdf5_folder <- tempdir()
settings <- GenoGAMSettings(hdf5Control = list(dir = hdf5_folder))

## C.
## register parallel backend (default is "the number of cores" - 2)
## Here we also use the SnowParam backend which is advised for larger data
## For yeast MulticoreParam should also do fine
BiocParallel::register(BiocParallel::SnowParam(worker = 2))

## D.
## specify chunk and overhang size. It is possible to skip this,
## but the automatic specification would prevent us from using
## multiple tiles in such a small example.
chunkSize <- 50000
overhangSize <- 1000

## E.
## the experiment design reflecting the underlying GAM
## x = position
design <- ~ s(x) + s(x, by = genotype)

## --------------------------------------------------------------------------
## build the GenoGAMDataSet with HDF5 backend
ggd <- GenoGAMDataSet(
  expDesign, directory = folder,
  chunkSize = chunkSize, overhangSize = overhangSize,
  design = design,
  settings = settings, hdf5 = TRUE
)

ggd

## --------------------------------------------------------------------------
## compute size factors
ggd <- computeSizeFactors(ggd)

ggd

## the data
assay(ggd)

## --------------------------------------------------------------------------
## compute model without parameter estimation to save time in vignette
result <- genogam(ggd, lambda = 4601, theta = 4.51)

result

## the fit and standard error
fits(result)
se(result)

## --------------------------------------------------------------------------
computeSignificance(result)

pval(result)

## --------------------------------------------------------------------------
idx <- 45000:55000
upper_input <- fits(result)$chrI[,"s(x)"][idx] + 2*se(result)$chrI[,"s(x)"][idx]
lower_input <- fits(result)$chrI[,"s(x)"][idx] - 2*se(result)$chrI[,"s(x)"][idx]

upper_ip <- fits(result)$chrI[,"s(x):genotype"][idx] + 2*se(result)$chrI[,"s(x):genotype"][idx]
lower_ip <- fits(result)$chrI[,"s(x):genotype"][idx] - 2*se(result)$chrI[,"s(x):genotype"][idx]

## ----plot, fig.cap = "The fit of input (upper) and IP (lower) for chrI:45000-55000", fig.width = 10, fig.height = 10----
par(mfrow = c(2,1))
plot(fits(result)$chrI[,"s(x)"][idx], type = "l", ylab = "log-fit input")
lines(upper_input, lty = 2)
lines(lower_input, lty = 2)

plot(fits(result)$chrI[,"s(x):genotype"][idx], type = "l", ylab = "log-fit IP")
lines(upper_ip, lty = 2)
lines(lower_ip, lty = 2)

## --------------------------------------------------------------------------
sessionInfo()

