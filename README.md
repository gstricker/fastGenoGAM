---
  title: "Modeling ChIP-Seq data with GenoGAM2: A Genome-wide generalized additive model"
  shorttitle: "fastGenoGAM"
  author: 
  - name: Georg Stricker
	affiliation: &id Technical University Munich
	email: georg.stricker@gmx.net
  - name: Julien Gagneur
	affiliation: *id
	email: gagneur@in.tum.de
  date: "`r format(Sys.Date(), '%m/%d/%Y')`"
  package: fastGenoGAM
  abstract: >
	Many genomic assays lead to noisy observations of a biological quantity of 
	interest varying along the genome. This is the case for ChIP-Seq, for which
	read counts reflect local protein occupancy of the ChIP-ed protein. 
	The fastGenoGAM package allows statistical analysis of genome-wide data with 
	smooth functions using generalized additive models. 
	<!-- It provides methods for -->
	<!-- the statistical analysis of ChIP-Seq data including inference of protein  -->
	<!-- occupancy, and pointwise and region-wise differential analysis as well as  -->
	<!-- peak calling with position-wise confidence bands. Estimation of dispersion  -->
	<!-- and smoothing parameters is performed by cross-validation. Scaling of  -->
	<!-- generalized additive model fitting to whole chromosomes is achieved by  -->
	<!-- parallelization over overlapping genomic intervals.  -->
	This vignette explains
	the use of the package for typical ChIP-Seq analysis workflow.
  output: 
    BiocStyle::html_document:
		toc_float: true
  bibliography: bibliog.bib
  vignette: >
    %\VignetteIndexEntry{Vignette Title}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}  
 ---

```{r setup, echo=FALSE, results="hide"}
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png",
                      message=FALSE, error=FALSE, warning=TRUE)
```

GenoGAM version: packageVersion("GenoGAM")

**Note:** if you use fastGenoGAM in published research, please cite:

> Stricker and Engelhardt, et al. (2017)
> GenoGAM: Genome-wide generalized additive models for ChIP-seq analysis
> *Bioinformatics*, **33**:15.
> [10.1093/bioinformatics/btx150](https://doi.org/10.1093/bioinformatics/btx150)

# TL;DR (with HDF5)

This is the brief version of the usual workflow of fastGenoGAM. It involves:

- Reading in data through `Rfunction{"GenoGAMDataSet"}` to get a `Robject{"GenoGAMDataSet"}` object. This is done with the HDF5 backend.
  
- Computing size factors with `Rfunction{"computeSizeFactors"}`
  
- Compute the model with `Rfunction{"genogam"}` to get the result `Robject{"GenoGAM"}` object

- compute position-wise log p-values

The dataset used is restricted to the first 100kb of the first yeast chromosome

First we set some need variables
```{r}
library(fastGenoGAM)

## specify folder and experiment design path
wd <- system.file("extdata/Set1", package='fastGenoGAM')
folder <- file.path(wd, "bam")
expDesign <- file.path(wd, "experimentDesign.txt")

## set HDF5 folder where the data will be stored
## Note, don't usually use /tmp because all your data and
## results get deleted otherwise later
hdf5_folder <- tempdir()
settings <- GenoGAMSettings(hdf5Control = list(dir = hdf5_folder))

## register parallel backend (default is "the number of cores" - 2)
## Here we also use the SnowParam backend which is advised for larger data
## For yeast MulticoreParam should also do fine
BiocParallel::register(BiocParallel::SnowParam(worker = 2))
```

Then we specify some optional parameters

```{r}
## specify chunk and overhang size. It is possible to skip this,
## but the automatic specification would prevent us from using
## multiple tiles in such a small example.
chunkSize <- 50000
overhangSize <- 1000

## the experiment design reflecting the underlying GAM
## x = position
design <- ~ s(x) + s(x, by = genotype)
```

And finally read in the data to obtain a `Robject{"GenoGAMDataSet"}` object.
Warnings about out-of-bound ranges can be ignored, as they occur during the 
shifting process when the genome gets tiled. It is trimmed accordingly.

```{r}
## build the GenoGAMDataSet with HDF5 backend
ggd <- GenoGAMDataSet(
  expDesign, directory = folder,
  chunkSize = chunkSize, overhangSize = overhangSize,
  design = design,
  settings = settings, hdf5 = TRUE
)

ggd
```

Compute Size factors

```{r}
## compute size factors
ggd <- computeSizeFactors(ggd)

ggd

## the data
assay(ggd)
```

Compute the model

```{r}
## compute model without parameter estimation to save time in vignette
result <- genogam(ggd, lambda = 4601, theta = 4.51)

result

## the fit and standard error
fits(result)
se(result)
```

Compute log p-values
```{r}
computeSignificance(result)

pval(result)
````

Plot results of the center 10kb, where both tiles are joined together.
```{r}
idx <- 45000:55000
upper_input <- fits(result)$chrI[,"s(x)"][idx] + 2*se(result)$chrI[,"s(x)"][idx]
lower_input <- fits(result)$chrI[,"s(x)"][idx] - 2*se(result)$chrI[,"s(x)"][idx]

upper_ip <- fits(result)$chrI[,"s(x):genotype"][idx] + 2*se(result)$chrI[,"s(x):genotype"][idx]
lower_ip <- fits(result)$chrI[,"s(x):genotype"][idx] - 2*se(result)$chrI[,"s(x):genotype"][idx]

par(mfrow = c(2,1))
plot(fits(result)$chrI[,"s(x)"][idx], type = "l")
lines(upper_input, lty = 2)
lines(lower_input, lty = 2)

plot(fits(result)$chrI[,"s(x):genotype"][idx], type = "l")
lines(upper_ip, lty = 2)
lines(lower_ip, lty = 2)
```

# FAQ

1. An error occurs during model fitting along those lines:

> error: matrix multiplication: incompatible matrix dimensions: 22333830147200x5360120267008000 and 4294972496x1

**Solution:** First, make sure you have all Armadillo dependencies installed correctly. See ![here](http://arma.sourceforge.net/download.html)

Second, the error is most likely related to the fact, that Armadillo is using 32bit matrices, thus causing problems for large matrices fastGenoGAM is using. The solution is to enable `ARMA_64BIT_WORD`, which is not enabled in RcppArmadillo by default. This should have been done during compilation, but if it fails for some reason you can do it manually with `#define ARMA_64BIT_WORD 1` in `my_R_Directory/lib/R/library/RcppArmadillo/include/RcppArmadilloConfig.h`

# Acknowledgments

We thank Alexander Engelhardt, Mathilde Galinier, Simon Wood, Herv\'e Pag\`es, and Martin Morgan for input in the development of fastGenoGAM

# Session Info

```{r}
sessionInfo()
```

# References

