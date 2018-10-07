#' Call peaks on a GenoGAM object
#'
#' Call narrow or broad peaks on the GenoGAM fit and computing significance, respectively 
#'
#' @param fit A \code{GenoGAM} object
#' @param smooth The name of the smooth, i.e. the colnames of the fitted object. By default
#' all are taken.
#' @param range A \code{GRanges} object specifying a range. By default the complete fit is
#' taken.
#' @param peakType The type of the peak (narrow or broad). Default is narrow, see details.
#' @param threshold The significance threshold. Keep in mind that the treshold depends on
#' the thresholdType. By default this is 0.05 for 'pvalue' and 0.1 for 'fdr'.
#' @param thresholdType The threshold type. Either 'fdr'(default) or 'pvalue'. If the
#' threshold is not provided it, will be set accordingly to the thresholdType.
#' @param maxgap For broad peaks only. The maximum gap between two broad peaks, that can be
#' tolerated in order to identify both as part of one broad peak. All broad peaks with
#' distances smaller or equal to the maxgap will be merged.
#' @param cutoff A seperate threshold for broad peaks. Since pointwise pvalues are
#' available. This threshold is used to identify all significantly high positions, which
#' then make up a broad peak.
#' @param minregion For broad peaks only. The minimum length of a broad peak. By default 1,
#' thus all found peaks are returned.
#' @return A list of data.tables of identified peaks. The different columns loosely resemble
#' the narrow and broad peak format (with different column order), such that it is easy to
#' write them to a 'narrowPeak', 'broadPeak' file (see \code{writeToBEDFile}).
#' @details Note, that broad peaks don't provide a specific highest location, but a region.
#' Whereas narrow peaks provide both. However, the borders of narrow peaks are not
#' necessarily informative but taken as +- 100bp around the peak summit. A function for a
#' more statistical estimation of the borders is being implemented. Also narrow peaks
#' provide an occupancy estimate at the peak position, while broad peaks give the average
#' occupancy accross the region.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
callPeaks <- function(fit, smooth = NULL, range = NULL,
                      peakType = c("narrow", "broad"), threshold = NULL,
                      thresholdType = c("fdr","pvalue"), maxgap = 500,
                      cutoff = 0.05, minregion = 1) {

    ## check choice parameters
    peakType <- match.arg(peakType)
    thresholdType <- match.arg(thresholdType)

    futile.logger::flog.info(paste0("Calling ", peakType, " peaks"))

    ## determine what type of object we are dealing with
    is_hdf5 <- is.HDF5(fit)
    if(class(fit) == "GenoGAMList") {
        is_split <- TRUE
    }
    else {
        if(class(fit) == "GenoGAM") {
            is_split <- FALSE
        }

        else {
            stop("Wrong class submitted")
        }
    }

    ## set smooth to all if not specified
    if(is.null(smooth)) {
        smooth <- colnames(fit)
    }
    
    ## set range to complete set if not specified
    if(is.null(range)) {
        range <- .extractGR(rowRanges(fit))
    }

    ## set threshold if not set
    if(is.null(threshold)) {
        threshold <- switch(thresholdType,
                            fdr = 0.1,
                            pvalue = 0.05)
    }

    ## compute background parameters
    background <- .computeBackground(fits(fit), se(fit), smooth)

    ## define grid over which to perform peak calling in parallel
    nranges <- 1:length(range)
    range_grid <- expand.grid(nranges, smooth)

    ## initialize result object
    peaks <- vector("list", length(smooth))
    names(peaks) <- smooth

    ## call peaks
    if(peakType == "narrow") {
        peakList <- BiocParallel::bplapply(1:nrow(range_grid), .callNarrowPeaks,
                                         grid = range_grid, fits = fits(fit),
                                         se = se(fit), rowRanges = rowRanges(fit),
                                         range = range, smooth = smooth,
                                         background = background)

        ## concatenate by smooths
        for(lev in unique(range_grid[,2])) {
            id <- which(range_grid[,2] == lev)
            peaks[[lev]] <- rbindlist(peakList[id]) 
        }
    }

    if(peakType == "broad") {
        peakList <- BiocParallel::bplapply(1:nrow(range_grid), .callBroadPeaks,
                                           grid = range_grid, fits = fits(fit),
                                           se = se(fit), rowRanges = rowRanges(fit),
                                           range = range, smooth = smooth,
                                           background = background, maxgap = maxgap,
                                           cutoff = cutoff)

        for(lev in unique(range_grid[,2])) {
            id <- which(range_grid[,2] == lev)
            peaks[[lev]] <- rbindlist(peakList[id])
            peaks[[lev]] <- peaks[[lev]][width >= minregion,]
        }
        
    }

    if(thresholdType == "pvalue") {
        for(lev in names(peaks)) {
            peaks[[lev]] <- peaks[[lev]][score >= -log(threshold),]
            peaks[[lev]] <- peaks[[lev]][order(score, decreasing = TRUE),]
        }
    }
    if(thresholdType == "fdr") {
        for(lev in names(peaks)) {
            peaks[[lev]] <- peaks[[lev]][fdr <= threshold,]
            peaks[[lev]] <- peaks[[lev]][order(fdr),]
        }
    }
    futile.logger::flog.info("DONE")
    return(peaks)
}

#' find local peaks in a signal
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.findPeaks <- function(x) {
    pks <- which(diff(sign(diff(x, na.pad=FALSE)),na.pad=FALSE) < 0) + 2
    return(pks)
}

#' find local minimas in a signal
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.findValleys <- function(x) {
    pks <- which(diff(sign(diff(x, na.pad=FALSE)),na.pad=FALSE) > 0) + 2
    return(pks)
}

## TODO: needs adjustment for hdf5 and split. This is for default only

#' compute parameters of a background distribution
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.computeBackground <- function(fit, se, smooth) {
    ## initializing the result matrix for smooths and mu0/var0 
    res <- data.frame(matrix(,length(smooth), 2), row.names = smooth)
    colnames(res) <- c("mu0", "var0")
        
    for(name in smooth) {
        all <- fit[,name]
        se_all <- se[,name]
        
        ## if genefilter package not present use median
        if(!requireNamespace("genefilter", quietly = TRUE)) {
            futile.logger::flog.info("genefilter package not found, median will be used instead of shorth.")
            mu0 <- median(all, na.rm=TRUE)
        }
        else {
            mu0 <- genefilter::shorth(all, na.rm=TRUE)
        }

        mu0_msg <- paste0("mu0 for ", name, " was computed as ", mu0)
        futile.logger::flog.debug(mu0_msg)
        
        left <- all[all <= mu0]
        right <- abs(left - mu0) + mu0
        new_data <- c(left, right)
        var0 <- mad(new_data, na.rm = TRUE)^2

        ## assign afterwards for readability
        res[name, 'mu0'] <- mu0
        res[name, 'var0'] <- var0

        var0_msg <- paste0("var0 for ", name, " was computed as ", var0)
        futile.logger::flog.debug(var0_msg)
    }
    return(res)
}

## ordered by zscore
.callNarrowPeaks <- function(iter, grid, fits, se, rowRanges, range, smooth, background) {
    r <- range[grid[iter, 1]] ## the range
    sx <- smooth[grid[iter, 2]]
    futile.logger::flog.debug(paste0("Calling peaks in region ", as.character(r)))

    mu0 <- background[sx, 'mu0']
    var0 <- background[sx, 'var0']

    ## find the region indices 
    idx <- S4Vectors::queryHits(IRanges::findOverlaps(rowRanges, r))
    startPos <- idx[1]

    ## find peaks
    region <- fits[idx, sx]
    peaks <- .findPeaks(region) + startPos - 1
    zscore <- (fits[peaks, sx] - mu0)/(sqrt(se[peaks, sx]^2 + var0))
    pvals <- -pnorm(-zscore, log.p = TRUE)
    dfpeaks <- data.table::data.table(chromosome = as.character(seqnames(r)),
                                      pos = peaks,
                                      zscore = zscore, score = pvals)
    
    ## for valleys
    valleys <- .findValleys(region) + startPos - 1
    vzscore <- (fits[valleys, sx] - mu0)/(sqrt(se[valleys, sx]^2 + var0))
    vzscore <- -vzscore
    dfvalleys <- data.table::data.table(pos = valleys, zscore = vzscore)

    ## compute FDR
    dfpeaks <- dfpeaks[order(zscore, decreasing = TRUE),]
    dfvalleys <- dfvalleys[order(zscore, decreasing = TRUE)]
    fdr <- sapply(dfpeaks$zscore, function(y) {
        sum(dfvalleys$zscore >= y)/sum(dfpeaks$zscore >= y)
    })
    ## adjust fdr for odd length of both vectors -->
    ## Number of peaks and valleys might differ by one, causing FDR become
    ## slightly greater than 1
    fdr <- fdr*(nrow(dfpeaks)/nrow(dfvalleys))
    dfpeaks$fdr <- fdr
    dfpeaks$summit <- exp(fits[peaks, sx])

    return(dfpeaks)
}

.callBroadPeaks <- function(iter, grid, fits, se, rowRanges, range, smooth, background,
                           maxgap, cutoff){
    r <- range[grid[iter, 1]] ## the range
    sx <- smooth[grid[iter, 2]]
    futile.logger::flog.debug(paste0("Calling peaks in region ", as.character(r)))

    mu0 <- background[sx, 'mu0']
    var0 <- background[sx, 'var0']

    ## find the region indices 
    idx <- S4Vectors::queryHits(IRanges::findOverlaps(rowRanges, r))
    startPos <- idx[1]

    ## compute zscore for all positions
    zscore <- (fits[idx, sx] - mu0)/(sqrt(se[idx, sx]^2 + var0))
    pvals <- -pnorm(-zscore, log.p = TRUE)
    pos <- start(r):end(r)

    signifPos <- pos[pvals >= -log(cutoff)]
    if(length(signifPos) == 0) {
      dt <- data.table::data.table()
    }
    else {
      diffs <- abs(diff(signifPos))
      breaks <- which(diffs > maxgap)
    
      starts <- c(signifPos[1], signifPos[breaks + 1])
      ends <- c(signifPos[breaks], signifPos[length(signifPos)])
      chrom <- seqnames(r)
      gr <- GenomicRanges::GRanges(chrom, IRanges(starts, ends))
      gp <- GenomicRanges::GPos(gr)

      ## non-significant position regions  < maxgap are incorporated into the
      ## broad peak. Thus we have more positions than signifPos and have to match
      ## and normalize them to start with 1 in order to use them as an index.
      indx <- match(pos(gp), pos) - pos[1] + 1
      gp$zscore <- zscore[indx]
      gp$pval <- pvals[indx]
      gp$estimate <- exp(fits[indx, sx])
      gp$region <- subjectHits(findOverlaps(gp, gr))
      
      ## compute significance for regions just like
      ## in differential binding: 1. Hochberg for a region-wise pvalue
      ## 2. BH for FDR of the regions between each other. 
      res <- data.table::data.table(as.data.frame(gr))
      dt <- data.table::data.table(as.data.frame(gp))
      pv <- dt[, min(p.adjust(exp(-pval), method="hochberg")), by = region]
      fpb <- dt[, mean(estimate), by = region]
      res$score[pv$region] <- pv$V1
      res$meanSignal[fpb$region] <- fpb$V1
      res$fdr = p.adjust(res$score, method="BH")
      res$score <- -log(res$score)
    }

    return(res)
}

#' Write peaks to BED6+3/4 format
#'
#' A function to write the data.table of peaks into a narrowPeaks or broadPeaks file
#'
#' @param peaks A data.table or data.frame of peaks as produced by 'callPeaks'
#' @param file A file name without suffix. It will be determined automatically. If no
#' file is given, it will be written to a generic 'peaks_[timestamp]' file in the
#' current working directory
#' @return Nothing. A narrowPeaks or broadPeaks file written to 'file'
#' @details Note, the narrow peak calling process does not yet implement any functionality
#' for estimating the start and end of a peak region. Thus the start and end is taken as -100
#' and +100 around the peak summit. This is mostly an arbitrary choice. A more statistical
#' approach is in development.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
writeToBEDFile <- function(peaks, file = NULL){
  writeBroadPeaks <- FALSE
  if("width" %in% names(peaks)) {
    writeBroadPeaks <- TRUE
  }
  if(is.null(file)) {
    timestamp <- gsub("-", "",strsplit(as.character(Sys.time()), split=" ")[[1]][1])
    file <- paste0("peaks_", timestamp)
  }
  if(writeBroadPeaks) {
    file <- paste0(file, ".broadPeak")
    writeToBroadPeaks(peaks, file)
  }
  else {
    file <- paste0(file, ".narrowPeak")
    writeToNarrowPeaks(peaks, file)
  }
  futile.logger::flog.info(paste0("File written to ", file.path(getwd(),file)))
}

#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
writeToNarrowPeaks <- function(peaks, file){
    ## so far start and stop are set arbitrary as 100bp around peak
    bp <- 100
    res <- peaks[,list(chrom = chromosome)]
    res$chromStart <- format(round(peaks$pos - bp), scientific = FALSE)
    res$chromEnd <- format(round(peaks$pos + bp), scientific = FALSE)
    res$name <- "."
    res$score <- "0"
    res$strand <- "."
    res$signalValue <- format(peaks$summit, scientific = FALSE)
    res$pValue <- format(-log10(exp(-peaks$pvalue)), scientific = FALSE)
    res$qValue <- format(-log10(peaks$fdr), scientific = FALSE)
    res$peak <- format(bp, scientific = FALSE)
    write.table(res, file = file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}

#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
writeToBroadPeaks <- function(peaks, file){
  res <- peaks[,list(chrom = seqnames, chromStart = start, chromEnd = end)]
  res$chromStart <- format(round(res$chromStart), scientific = FALSE)
  res$chromEnd <- format(round(res$chromEnd), scientific = FALSE)
  res$name <- "."
  res$score <- "0"
  res$strand <- "."
  res$signalValue <- format(peaks$meanSignal, scientific = FALSE)
  res$pValue <- format(-log10(exp(-peaks$score)), scientific = FALSE)
  res$qValue <- format(peaks$fdr, scientific = FALSE)
  write.table(res, file = file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}
