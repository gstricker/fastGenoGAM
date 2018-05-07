###############################
## p-value computation
###############################

#' computeSignificance
#'
#' The function computes positionwise p-values
#' 
#' @param gg A GenoGAM object.
#' @param log.p Should values be returned in log scale?
#' @param nquantiles The number of quantiles for base level computation.
#' Only effective if the GenoGAM object is of class GenoGAMList or
#' the underlying data is stored in HDF5 format. See details.
#' @return An updated GenoGAM object, where the p-value slot is added.
#' @details In case of GenoGAMList the base level is computed as the median over
#' a list of chromosome data. This is achieved by computing a number of quantiles
#' for each chromosome and then compute the median on the complete set of quantiles
#' over all chromosomes. Similarly the same is done for the underlying HDF5 data,
#' but possibly with an even larger set of quantiles. That is because the data is read and
#' quantiles computed for each block to keep memory footprint low. Unless the organisms is
#' small, blocks have smaller size than chromosomes. Note, that this is an approximative
#' approach and the median is not necessarily closer to the overall 'real' median with
#' larger number of quantiles (unless in the extreme case, where each quantile is one data point).
#' 5 quantiles (the default setting) is a number that gave good results in practice.
#'
#' Note, that in case the data is stored in HDF5 format, the pvalue 'group' is added on hard drive.
#' That is, unlike any other function in R, where the input object is not changed, it actually
#' is in this case. If one wishes to have HDF5 data without the pvalue 'group', one has to
#' backup the HDF5 files prior to computation.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
computeSignificance <- function(gg, log.p = TRUE, nquantiles = 5) {

    if(log.p) {
        futile.logger::flog.info("Computing positionwise log p-values")
    }
    else {
        futile.logger::flog.info("Computing positionwise p-values")
    }

    .pvals(gg, log.p, nquantiles)
    
    futile.logger::flog.info("DONE")
}

.pvals <- function(gg, log.p = FALSE, nquantiles = 5) {
    hdf5 <- is.HDF5(gg)
    split <- (class(gg) == "GenoGAMList")

    if(split){
        
        if(hdf5) {
            futile.logger::flog.debug("Data is in HDF5 format and split by chromosome")
            .pvals_hdf5_split(gg, log.p, nquantiles)
        }
        else {
            futile.logger::flog.debug("Data is split by chromosome")
            .pvals_split(gg, log.p, nquantiles)
        }
    }
    else {
        if(hdf5) {
            futile.logger::flog.debug("Data is in HDF5 format")
            .pvals_hdf5(gg, log.p, nquantiles)
        }
        else {
            .pvals_default(gg, log.p)
        }
    }
    invisible(NULL)
}

## compute pvalue for GenoGAM object without HDF5 and no split
.pvals_default <- function(gg, log.p = FALSE) {
    tracks <- colnames(gg)
    pvals <- lapply(tracks, function(id) {
        ## base <- median(fits(gg)[[id]]) ## get the logarithmic base level of the background
        base <- 0
        ## futile.logger::flog.debug(paste("The base level for track", id, "computed as:", base))
        
        res <- 2*pnorm(base, mean = abs(fits(gg)[[id]]), sd = se(gg)[[id]], log.p = log.p)
        return(res)
    })
    df <- DataFrame(data.frame(pvals))
    names(df) <- tracks

    ## add pvals to assay fields. SummarizedExperiment has some control functions
    ## guarding the structure of the data, so we have to directly plug it in to
    ## avoid copying
    gg@assays$data@listData$pval <- df
    ## return(gg)
}

## compute pvalue for GenoGAM object with HDF5 and no split
.pvals_hdf5 <- function(gg, log.p = FALSE, nquantiles = 5) {
    tracks <- colnames(gg)

    ## make HDF5 group in existing HDF5 file for p-values
    dims <- dim(assay(gg))
    
    seedFile <- assay(gg)@seed@seed@filepath
    h5file <- rhdf5::H5Fopen(seedFile)
    chunks <- getHDF5DumpChunkDim(dims, "double")
    h5name <- "pvals"
    tryCatch(
        .createH5Dataset(h5file, name = h5name, settings = gg@settings, d = dims,
                         chunk = chunks, type = "H5T_IEEE_F32LE"),
        error = function(e) {
            warning("'pvals' group already exists. Overwriting.")
        })
    ## close file to save initial zero matrix
    rhdf5::H5Fclose(h5file)
    
    ## reopen to write again
    h5file <- rhdf5::H5Fopen(seedFile)
    pvals <- lapply(1:length(tracks), function(id) {
        sp_fits <- fits(gg)[,id]
        sp_ses <- se(gg)[,id]
        ## base <- .compute_base(sp_fits, nquantiles)
        base <- 0
        futile.logger::flog.debug(paste("The base level for track", tracks[id], "computed as:", base))
        
        .hdf5_block_pval(x = sp_fits, se = sp_ses, base = base,
                         h5file = h5file, name = h5name, id = id,
                         log.p = log.p, nquantiles = nquantiles)
    })

    ## add pvals to assay fields. SummarizedExperiment has some control functions
    ## guarding the structure of the data, so we have to directly plug it in to
    ## avoid copying
    h5 <- HDF5Array::HDF5Array(seedFile, h5name)
    gg@assays$data@listData$pval <- h5
    return(gg)
}

## compute pvalue for GenoGAM object with HDF5 and no split
.pvals_hdf5_split <- function(gg, log.p = FALSE, nquantiles = 5) {
    tracks <- colnames(gg)
    sp_fits <- fits(gg)
    sp_ses <- se(gg)
    chroms <- names(sp_fits)

    ## compute base level
    ## base <- .compute_base(sp_fits, nquantiles)
    base <- 0
    if(futile.logger::flog.threshold() == "DEBUG") {
        for(ii in 1:length(base)) {
            futile.logger::flog.debug(paste("The base level for track", tracks[ii], "computed as:", base[[ii]]))
        }
    }

    
    lapply(chroms, function(y) {
        futile.logger::flog.debug(paste("Computing p-values for chromosome", y))
        ## make HDF5 group in existing HDF5 file for p-values
        dims <- dim(assay(gg@data[[y]]))
    
        seedFile <- assay(gg@data[[y]])@seed@seed@filepath
        h5file <- rhdf5::H5Fopen(seedFile)
        chunks <- assay(gg@data[[y]])@seed@chunkdim
        h5name <- "pvals"
        tryCatch(
            .createH5Dataset(h5file, name = h5name, settings = gg@settings, d = dims,
                             chunk = chunks, type = "H5T_IEEE_F32LE"),
            error = function(e) {
                warning("'pvals' group already exists. Overwriting.")
            })
        ## close file to save initial zero matrix
        rhdf5::H5Fclose(h5file)
        
        ## reopen to write again
        h5file <- rhdf5::H5Fopen(seedFile)
        pvals <- lapply(1:length(tracks), function(id) { 
            .hdf5_block_pval(x = sp_fits[[y]][,id], se = sp_ses[[y]][,id],
                             base = base[[id]], h5file = h5file,
                             name = h5name, id = id, log.p = log.p,
                             nquantiles = nquantiles)
        })

        ## add pvals to assay fields. SummarizedExperiment has some control functions
        ## guarding the structure of the data, so we have to directly plug it in to
        ## avoid copying
        h5 <- HDF5Array::HDF5Array(seedFile, h5name)
        gg@data[[y]]@assays$data@listData$pval <- h5
    })

    ## return(gg)
}

.median <- function(block, nquantiles = 5) {
    
    APPLY <- function(block, nquantiles) {
        tmp <- as.vector(block, mode = "numeric")
        res <- quantile(tmp, probs = seq(0 , 1, length.out = nquantiles))
        return(res)
    }

    COMBINE <- function(x) {
        median(x)
    }

    return(list(APPLY = APPLY, COMBINE = COMBINE))
}

## make equal size grid to avoid quantiles on small regions which are
## then over-represented in the complete array
.makeGrid <- function(x) {
    max_block_size <- DelayedArray:::get_max_block_length(type(x))
    xsize <- nrow(x)
    mult <- as.integer(ceiling(xsize/max_block_size))
    optimal_block_size <- as.integer(xsize/mult)
    grid <- DelayedArray::defaultGrid(x, optimal_block_size)
    
    return(grid)
}

.block_APPLY <- function(x, fun, nquantiles = 5) {

    ## if it's a list, then every element is most likely a matrix --> go by column
    if(is(x, "list")) {
        ntracks <- ncol(x[[1]])
        res <- lapply(1:ntracks, function(id) {
            qn <- sapply(x, function(y) {
                grid <- .makeGrid(y[,id])
                nblock <- length(grid)
                sqn <- numeric(nblock * nquantiles)
                for (b in seq_len(nblock)) {
                    ## take block and compute the 1% quantiles of it
                    block <- DelayedArray:::extract_block(y[,id], grid[[b]])
                    idx <- (1:nquantiles) + nquantiles * (b - 1)
                    sqn[idx] <- fun$APPLY(block, nquantiles)
                }
                return(sqn)
            })
            qn <- unlist(qn)
            res <- fun$COMBINE(qn)
            return(res)
        })
    }

    else {
        grid <- .makeGrid(x)
        nblock <- length(grid)

        qn <- numeric(nblock * nquantiles)
        for (b in seq_len(nblock)) {
            ## take block and compute the 1% quantiles of it
            block <- DelayedArray:::extract_block(x, grid[[b]])
            idx <- (1:nquantiles) + nquantiles * (b - 1)
            qn[idx] <- fun$APPLY(block, nquantiles)
        }
        res <- fun$COMBINE(qn)
    }
    
    return(res)
}

.compute_base <- function(x, nquantiles = 5) {
    ## compute median
    base <- .block_APPLY(x, .median(), nquantiles)

    return(base)
}

.hdf5_block_pval <- function(x, se, base, h5file, name, id, log.p = FALSE, nquantiles = 5) {

    grid <- .makeGrid(x)
    nblock <- length(grid)

    ## now compute pvalues
    for (b in seq_len(nblock)) {
        subgrid <- grid[[b]]
        xblock <- DelayedArray:::extract_block(x, subgrid)
        seblock <- DelayedArray:::extract_block(se, subgrid)
        xtmp <- as.vector(xblock, mode = "numeric")
        setmp <- as.vector(seblock, mode = "numeric")

        res <- 2*pnorm(base, mean = abs(xtmp), sd = setmp, log.p = log.p)
        idx <- ranges(subgrid)[1]
        rhdf5::h5write(as.matrix(res, ncol = 1), file = h5file, name = name,
                       index = list(start(idx):end(idx), id))
    }
    invisible(NULL)
}

## compute pvalue for GenoGAM object without HDF5, but with Split
.pvals_split <- function(gg, log.p = FALSE, nquantiles = 5) {
    tracks <- colnames(gg)
    chroms <- names(assays(gg))
    len <- length(chroms)

    ## first get median of medians
    ## chromosome dependent medians vary too much, especially for small organisms
    ## base <- lapply(tracks, function(id) {
    ##     sqn <- numeric(nquantiles * len)
    ##     for(ii in 1:len) {
    ##         idx <- (1:nquantiles) + nquantiles * (ii - 1)
    ##         sqn[idx] <- quantile(fits(gg)[[ii]][[id]], probs = seq(0 , 1, length.out = nquantiles))
    ##     }
    ##     res <- median(sqn)
    ##     futile.logger::flog.debug(paste("The base level for track", id, "computed as:", res))
    ##     return(res)
    ## })
    base <- list(0, 0)
    names(base) <- tracks

    ## now compute p-values
    lapply(chroms, function(y) {
        futile.logger::flog.debug(paste("Computing p-values for chromosome", y))
        pvals <- lapply(tracks, function(id) {
            res <- -2*pnorm(base[[id]], mean = abs(fits(gg)[[y]][[id]]), sd = se(gg)[[y]][[id]], log.p = log.p)
            return(res)
        })
        df <- DataFrame(data.frame(pvals))
        names(df) <- tracks

        ## add pvals to assay fields. SummarizedExperiment has some control functions
        ## guarding the structure of the data, so we have to directly plug it in to
        ## avoid copying
        gg@data[[y]]@assays$data@listData$pval <- df
        invisible(NULL)
    })
    ## return(gg)
}
