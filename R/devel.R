## library(fastGenoGAM)

## config <- "/s/project/coreProm/analysis/diffBinding/config.txt"
## folder <- "/s/project/coreProm/Michi/thornton/align_STAR"

## settings <- GenoGAMSettings(hdf5Control = list(dir = "/s/project/coreProm/hdf5/sacCer3"))
## ggd <- GenoGAMDataSet(config, chunkSize = 140000, overhangSize = 1000, design = ~ s(x) + s(x, by = genotype),
##                       directory = folder, fromHDF5 = TRUE, split = TRUE, settings = settings, ignoreM = TRUE)
## ggd <- computeSizeFactors(ggd)
## gg <- GenoGAM(ggd = ggd, family = "nb", fromHDF5 = TRUE, split = TRUE)

## ## pvalue computation
## .pvals <- function(gg, log.p = FALSE) {
##     hdf5 <- is.HDF5(gg)
##     split <- (class(gg) == "GenoGAMList")

##     if(split){
##         if(hdf5) {
##             lapply(names(assay(gg)), function(chrom) {

##                 dims <- dim(assay(gg)[[chrom]])

##                 ## make HDF5 group in existing HDF5 file for p-values
##                 seedFile <- assay(gg)[[1]]@seed@file
##                 h5file <- rhdf5::H5Fopen(seedFile)
##                 .createH5Dataset(h5file, name = "pvals", settings = settings, d = dims,
##                                  chunk = dims, type = "H5T_IEEE_F32LE")
##                 rhdf5::H5Fclose(h5file)

##                 sapply(1:dims[2], function(y) {
##                     fits <- assays(gg)[[chrom]][["fits"]][,y]
##                     ses <- assays(gg)[[chrom]][["se"]][,y]
##                     pval <- 2*pnorm(0, mean = fits, sd = ses, log.p = log.p) ## own Cpp function needed here
##                 })
##             })
##         }
##         else {
            
##         }
##     }
##     else {
##         if(hdf5) {
            
##         }
##         else {
            
##         }
##     }




    
    
##     cols <- names(gg@fits)
##     colindx <- which(unlist(regexec("se\\.", names(gg@fits))) < 0)

##     res <- data.frame(matrix(NA, nrow(gg@fits), length(colindx)))
##     for(ii in 1:length(colindx)) {
##         colname <- cols[colindx[ii]]
##         secolname <- paste("se", colname, sep = ".")
##         background <- (regexec(":", colname)[[1]] < 0)
##         if(background) {
##             intercept <- median(gg@fits[,colname], na.rm = TRUE)
##             fitMean <- abs(gg@fits[,colname] - intercept)
##         }
##         else fitMean <- abs(gg@fits[,colname])
##         res[,ii] <- 2*pnorm(0, mean = fitMean, sd = gg@fits[, secolname], log.p = log.p)
##     }
##     names(res) <- paste("pvalue", cols[colindx], sep = ".")
##     return(res)
## }

## .result <- function(ggd, log.p = FALSE) {
##     if(length(ggd) > 0) {
##         pvals <- .pvals(ggd, log.p)
##         slot(gg, "fits") <- cbind(slot(gg, "fits"), pvals)
##     }
##     return(gg)
## }

## #' Compute significance.
## #'
## #' Based on the model fits this functions computes pointwise pvalues.
## #'
## #' @param gg A fitted GenoGAM object.
## #' @param log.p Should pvalues be returned in log scale?
## #' @return A GenoGAM object which fits has been updated by the pvalue columns.
## #' @examples
## #' ggd <- makeTestGenoGAM()
## #' ggd <- computeSignificance(ggd)
## #' head(getFits(ggd))
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## #' @export
## computeSignificance <- function(gg, log.p = FALSE) {
##     .result(gg, log.p = log.p)
## }
