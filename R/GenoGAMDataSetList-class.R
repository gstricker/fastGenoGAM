## ====================
## GenoGAMDatasetList class
## ====================
#' @include GenoGAMSettings-class.R
NULL

#' GenoGAMDataSetList
#'
#' The GenoGAMDataSetList class contains the pre-processed raw data and
#' additional slots that define the input and framework for the model.
#' It extends upon the idea of the GenoGAMDataSet class to make it possible
#' to store genomes and data of size > 2^31 (maximum size of integers in R).
#' Thus the only difference to a GenoGAMDataSet is the arrangement as a
#' list of RangedSummarizedExperiments under the hood. On the surface is
#' does still behave like a GenoGAMDataSet. It is not intended to be used
#' by the user. For more information check the GenoGAMDataSet class documentation.
#' 
#' @slot settings The global and local settings that will be used to compute the
#' model.
#' @slot design The formula describing how to evaluate the data. See details.
#' @slot sizeFactors The normalized values for each sample. A named numeric vector.
#' @slot index A GRanges object representing an index of the ranges defined 
#' on the genome. Mostly used to store tiles.
#' @slot ggd A list of RangedSummarizedExperiment objects
#' @slot id A GRanges object keeping the identifiers assigning the regions to the
#' respective list elements
#' @name GenoGAMDataSetList-class
#' @rdname GenoGAMDataSetList-class
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
setClass("GenoGAMDataSetList",
         slots = list(settings = "GenoGAMSettings",
                      design = "formula", sizeFactors = "numeric",
                      index = "GRanges", ggd = "list", id = "GRanges"),
         prototype = list(settings = GenoGAMSettings(),
                          design = ~ s(x), sizeFactors = numeric(), 
                          index = GenomicRanges::GRanges(),
                          ggd = list(), id = GenomicRanges::GRanges()))

## Validity
## ========

.validateGGDType <- function(object) {
    if(class(object@ggd) != "list") {
        return("'ggd' must be a list object and all elements must be of class RangedSummarizedExperiment")
    }
    NULL
}

.validateIDType <- function(object) {
    if(class(object@id) != "GRanges") {
        return("'id' must be a GRanges object")
    }
    NULL
}

.validateGGDLChromosomes <- function(object) {
    cindex <- GenomeInfoDb::seqlevels(object@index)
    seqlev <- lapply(object@ggd, function(y) {
        GenomeInfoDb::seqlevels(SummarizedExperiment::rowRanges(y))
    })
    cobject <- unique(unlist(seqlev))
    if(!all(cindex %in% cobject)) {
        return("Different chromosomes for data and index objects.")
    }
    NULL
}

## general validate function
.validateGenoGAMDataSetList <- function(object) {
    c(.validateSettingsType(object), .validateDesignType(object),
      .validateSFType(object), .validateIndexType(object),
      .validateGGDLChromosomes(object),
      .validateGGDType(object),
      .validateIDType(object))
}

S4Vectors::setValidity2("GenoGAMDataSetList", .validateGenoGAMDataSetList)

## Constructor
## ===========

#' GenoGAMDataSetList constructor.
#'
#' GenoGAMDataSetList is the constructor function for the GenoGAMDataSetList-class.
#'
#' @aliases dim length seqlengths seqlevels colData rowRanges colnames getIndex
#' tileSettings dataRange getChromosomes getTileSize getChunkSize getOverhangSize
#' getTileNumber design design<- sizeFactors sizeFactors<- getChunkSize<-
#' getTileSize<- getOverhangSize<- getTileNumber<-
#' @param ... The slots and their respective values
#' @return An object of class GenoGAMDataSetList
#' @name GenoGAMDataSetList
#' @rdname GenoGAMDataSetList-class
GenoGAMDataSetList <- function(...) {
    return(new("GenoGAMDataSetList", ...))
}


#' Make an example /code{GenoGAMDataSet}
#'
#' @param sim Use simulated data (TRUE) or test data from a real experiment
#' @return A /code{GenoGAMDataSet} object
makeTestGenoGAMDataSetList <- function() {

    k <- 10000
    sinCurve <- sin(seq(-7, 5, length.out = k)) + 1
    ip <- rnbinom(k, size = 2, mu = sinCurve/max(sinCurve))
    sinCurve <- c(sin(seq(-7, -1, length.out = k/2)) + 1, runif(k/2, 0, 0.2))
    background <- rnbinom(k, size = 2, mu = sinCurve/max(sinCurve)/2)
    chroms <- c("chrX", "chrY", "chrZ")
    gr <- GenomicRanges::GPos(GenomicRanges::GRanges(chroms,
                                                     IRanges::IRanges(c(1,1,1), c(k, k, k))))
    GenomeInfoDb::seqlengths(gr) <- c(1e6, 2e6, 1.5e6)

    ## colData
    coldf <- S4Vectors::DataFrame(experiment = c(0, 1))
    rownames(coldf) <- c("input", "IP")
    
    selist <- lapply(chroms, function(chr) {
        add <- sample(0:3, 2)
        background <- background + add[1]
        ip <- ip + add[2]
        df <- S4Vectors::DataFrame(input = background, IP = ip)
        se <- SummarizedExperiment::SummarizedExperiment(rowRanges = gr[GenomeInfoDb::seqnames(gr) == chr,],
                                                         assays = list(df), colData = coldf)
    })
    names(selist) <- chroms
    id <- .extractGR(gr)
    id$id <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(gr)))

    ## make tiles
    l <- list(chromosomes = gr,
              chunkSize = 2000,
              overhangSize = 250)
    tiles <- .makeTiles(l)

    ## size factors
    sf <- rep(0, length(colnames(selist[[1]])))
    
    ggd <- GenoGAMDataSetList(ggd = selist, id = id, design = ~ s(x),
                              index = tiles, sizeFactors = sf)
    
    design(ggd) <- ~ s(x) + s(x, by = experiment)
        
    return(ggd)
}

## Check function
## ===============

#' @noRd
setMethod("checkObject", "GenoGAMDataSetList", function(object) {
    .checkGenoGAMDataSet(object)
})

## Accessors
## =========

#' @describeIn GenoGAMDataSetList Get the dimension of the object
setMethod("dim", "GenoGAMDataSetList", function(x) {
    dims <- sapply(ggd@ggd, dim)
    res <- c(sum(dims[1,]), max(dims[2,]))
    return(res)
})

#' @describeIn GenoGAMDataSetList The length of the object
setMethod("length", "GenoGAMDataSetList", function(x) {
    dim(x)[1]
})

#' @describeIn GenoGAMDataSetList The seqlengths of the object
setMethod("seqlengths", "GenoGAMDataSetList", function(x) {
    return(GenomeInfoDb::seqlengths(rowRanges(x)[[1]]))
})

#' @describeIn GenoGAMDataSetList The seqlevels of the object
setMethod("seqlevels", "GenoGAMDataSetList", function(x) {
    return(GenomeInfoDb::seqlevels(rowRanges(x)[[1]]))
})

##' @describeIn GenoGAMDataSetList get colData from the first element of the
##' SummarizedExperiment list
setMethod("colData", "GenoGAMDataSetList", function(x, ...) {
    return(colData(x@ggd[[1]]))
})

##' @describeIn GenoGAMDataSetList get a list of rowRanges from the
##' GenoGAMDataSetList object
setMethod("rowRanges", "GenoGAMDataSetList", function(x, ...) {
    lapply(x@ggd, rowRanges)
})

##' @describeIn GenoGAMDataSetList get a list of assays from the
##' GenoGAMDataSetList object
setMethod("assay", "GenoGAMDataSetList", function(x, ...) {
    lapply(x@ggd, assay)
})

##' @describeIn GenoGAMDataSetList get colnames from the first element of the
##' SummarizedExperiment list
setMethod("colnames", "GenoGAMDataSetList", function(x) {
    return(colnames(x@ggd[[1]]))
})

##' @describeIn GenoGAMDataSetList accessor to the index slot
setMethod("getIndex", "GenoGAMDataSetList", function(object) {
    return(slot(object, "index"))
})

##' @describeIn GenoGAMDataSetList The accessor to the list of settings, that
##' were used to generate the tiles.
setMethod("tileSettings", "GenoGAMDataSetList", function(object) {
    S4Vectors::metadata(getIndex(object))
})

#' @noRd
.rowRangesFromList <- function(object) {
    res <- lapply(object, function(x) {
        .extractGR(rowRanges(x))
    })
    return(do.call("c", unname(res)))
}

##' @describeIn GenoGAMDataSetList The actual underlying GRanges showing the range of the data.
setMethod("dataRange", "GenoGAMDataSetList", function(object) {
    res <- .rowRangesFromList(object@ggd)
})

##' @describeIn GenoGAMDataSetList A GRanges object representing the chromosomes
##' or chromosome regions on which the model will be computed
setMethod("getChromosomes", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$chromosomes
})

##' @describeIn GenoGAMDataSetList The size of the tiles
setMethod("getTileSize", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$tileSize
})

##' @describeIn GenoGAMDataSetList The size of the chunks
setMethod("getChunkSize", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$chunkSize
})

##' @describeIn GenoGAMDataSetList The size of the overhang (on one side)
setMethod("getOverhangSize", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$overhangSize
})

##' @describeIn GenoGAMDataSetList The total number of tiles
setMethod("getTileNumber", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$numTiles
})

##' @describeIn GenoGAMDataSetList Access to the design slot.
setMethod("design", "GenoGAMDataSetList", function(object) {
    slot(object, "design")
})

##' @describeIn GenoGAMDataSetList Replace method of the design slot.
setReplaceMethod("design", "GenoGAMDataSetList", function(object, value) {
    newCols <- as.vector(na.omit(.getVars(value)))
    if(!all(newCols %in% colnames(colData(object)))) {
        futile.logger::flog.error("'by' variables could not be found in colData")
        stop("'by' variables could not be found in colData")
    }
    slot(object, "design") <- value
    return(object)
})

##' @describeIn GenoGAMDataSetList Access to the sizeFactors slot
setMethod("sizeFactors", "GenoGAMDataSetList", function(object) {
    sf <- slot(object, "sizeFactors")
    names(sf) <- colnames(object)
    return(sf)
})

##' @describeIn GenoGAMDataSetList Replace method of the sizeFactors slot
setReplaceMethod("sizeFactors", "GenoGAMDataSetList", function(object, value) {
    slot(object, "sizeFactors") <- value
    return(object)
})


## change Settings
## ===============

##' @describeIn GenoGAMDataSetList Replace method of the chunkSize parameter,
##' that triggers a new computation of the tiles based on the new chunk size.
##' @noRd
setReplaceMethod("getChunkSize", signature = c("GenoGAMDataSetList", "numeric"),
                 function(object, value) {
                     settings <- tileSettings(object)
                     settings$chunkSize <- value
                     settings$tileSize <- value + 2*settings$overhangSize
                     newIndex <- .makeTiles(settings)
                     slot(object, "index") <- newIndex
                     return(object)
                 })

##' @describeIn GenoGAMDataSetList Replace method of the tileSize parameter,
##' that triggers a new computation of the tiles based on the new tile size.
##' @noRd
setReplaceMethod("getTileSize", signature = c("GenoGAMDataSetList", "numeric"),
                 function(object, value) {
                     settings <- tileSettings(object)
                     settings$tileSize <- value
                     settings$chunkSize <- value - 2*settings$overhangSize
                     newIndex <- .makeTiles(settings)
                     slot(object, "index") <- newIndex
                     return(object)
                 })

##' @describeIn GenoGAMDataSetList Replace method of the overhangSize parameter,
##' that triggers a new computation of the tiles based on the new overhang size.
##' @noRd
setReplaceMethod("getOverhangSize", signature = c("GenoGAMDataSetList", "numeric"),
                 function(object, value) {
                     settings <- tileSettings(object)
                     settings$overhangSize <- value
                     settings$tileSize <- settings$chunkSize + 2*value
                     newIndex <- .makeTiles(settings)
                     slot(object, "index") <- newIndex
                     return(object)
                 })

##' @describeIn GenoGAMDataSetList Replace method of the tileNumber parameter,
##' that triggers a new computation of the tiles based on the new number of tiles.
##' @noRd
setReplaceMethod("getTileNumber", signature = c("GenoGAMDataSetList", "numeric"),
                 function(object, value) {
                     settings <- tileSettings(object)
                     size <- settings$chunkSize*settings$numTiles
                     if(size > sum(width(dataRange(object)))) {
                         warning("The settings indicated a longer total genome size than actually present. it was trimmed accordingly.")
                         size <- sum(width(dataRange(object)))
                     }
                     settings$chunkSize <- round(size/value)
                     settings$tileSize <- value + 2*settings$overhangSize
                     newIndex <- .makeTiles(settings)
                     slot(object, "index") <- newIndex
                     return(object)
                 })

## Subsetting
## ==========
#' Subset method for GenoGAMDataSetList
#'
#' @details
#' Those are various methods to subset the GenoGAMDataSetList object.
#' By logical statement or GRanges overlap. The '[' subsetter is
#' just a short version of 'subsetByOverlaps'. The double brackets '[['
#' offer a subset based on tiles.
#'
#' @aliases subsetByOverlaps '[' '[['
#' @param x,query A GenoGAMDataSet object.
#' @param subject,i A GRanges object. In case of subsetting by double brackets
#' 'i' is the index of the tile.
#' @param maxgap,minoverlap Intervals with a separation of 'maxgap' or
#' less and a minimum of 'minoverlap' overlapping positions, allowing for
#' 'maxgap', are considered to be overlapping. 'maxgap' should
#' be a scalar, non-negative, integer. 'minoverlap' should be a scalar,
#' positive integer.
#' @param type By default, any overlap is accepted. By specifying the 'type'
#' parameter, one can select for specific types of overlap. The types correspond
#' to operations in Allen's Interval Algebra (see references in). If \code{type}
#' is \code{start} or \code{end}, the intervals are required to have matching
#' starts or ends, respectively. While this operation seems trivial, the naive
#' implementation using \code{outer} would be much less efficient. Specifying
#' \code{equal} as the type returns the intersection of the \code{start} and
#' \code{end} matches. If \code{type} is \code{within}, the query interval must
#' be wholly contained within the subject interval. Note that all matches must
#' additionally satisfy the \code{minoverlap} constraint described above.
#'
#' The \code{maxgap} parameter has special meaning with the special
#' overlap types. For \code{start}, \code{end}, and \code{equal}, it specifies
#' the maximum difference in the starts, ends or both, respectively. For
#' \code{within}, it is the maximum amount by which the query may be wider
#' than the subject.
#' @param invert If TRUE, keep only the query ranges that do _not_ overlap
#' the subject.
#' @param ... Further arguments. Mostly a logical statement
#' in case of the 'subset' function. Note that the columnnames
#' for chromosomes and positions are: 'seqnames' and 'pos'.
#' @references
#' Allen's Interval Algebra: James F. Allen: Maintaining knowledge
#' about temporal intervals. In: Communications of the ACM.
#' 26/11/1983. ACM Press. S. 832-843, ISSN 0001-0782
#' @return A subsetted GenoGAMDataSetList object.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @rdname GenoGAMDataSetList-subsetting
setMethod("subset", "GenoGAMDataSetList", function(x, ...) {
    if(all(dim(x) == c(0, 0))) return(x)
    
    settings <- slot(x, "settings")
    design <- design(x)
    sf <- sizeFactors(x)

    ## make ggd slot
    ## subset all SummarizedExperiments
    se <- lapply(x@ggd, function(se) {
        res <- subset(se, ...)
        if(length(res) == 0) {
            res <- NULL
        }
        return(res)
    })

    ## reduce the list to !NULL
    se <- se[!vapply(se, is.null, logical(1))]

    ## Check dimensions
    dims <- sapply(se, dim)
    dims <- c(sum(dims[1,]), max(dims[2,]))

    ## make index slot
    ## make empty index if data is empty otherwise trim accordingly
    if(any(dims == 0)) {
        index <- GenomicRanges::GRanges()
    }
    else {
        se <- lapply(se, function(y) {
            ## set correct seqinfo
            GenomeInfoDb::seqlevels(rowRanges(y), pruning.mode="coarse") <- GenomeInfoDb::seqlevelsInUse(rowRanges(y))
            return(y)
        })
        slot(settings, "chromosomeList") <- GenomeInfoDb::seqlevels(se[[1]])
        index <- .subsetIndex(se, getIndex(x))
    }

    ## make id slot
    splitid <- .rowRangesFromList(se)
    splitid$id <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(splitid)))
    
    ggd <- new("GenoGAMDataSetList", settings = settings,
               design = design, sizeFactors = sf, index = index,
               ggd = se, id = splitid)
    return(ggd)
})

.subsetIndexGDDL <- function(se, index) {
    ## get the data range
    gpCoords <- .rowRangesFromList(se)
   
    l <- S4Vectors::metadata(index)
    l$chromosomes <- gpCoords
    minWidth <- min(width(gpCoords))
    if(minWidth < l$tileSize) {
        l$tileSize <- minWidth
        l$chunkSize <- minWidth - 2*l$overhangSize
    }
    indx <- .makeTiles(l)
    GenomeInfoDb::seqinfo(indx) <- GenomeInfoDb::seqinfo(gpCoords)
    
    return(indx)
}

#' underlying function to subset by overlaps
.subsetByOverlapsGDDL <- function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"),
                   invert = FALSE, ...) {

    ## need to make even widths, otherwise subsetByOverlaps does it with a warning
    if(any((width(subject) %% 2) == 1)) {
        futile.logger::flog.info("Some subset ranges have odd widths. Rounding to the next even number.")
        idx <- which((width(subject) %% 2) == 1)
        width(subject)[idx] <- width(subject)[idx] + 1
    }          
    settings <- slot(query, "settings")
    design <- design(query)
    sf <- sizeFactors(query)

    ## iterate over all SummarizedExperiments and subset
    se <- lapply(query@ggd, function(se) {
        res <- subsetByOverlaps(se, subject, maxgap = maxgap,
                                minoverlap = minoverlap,
                                type=type, invert = invert, ...)
        if(length(res) == 0) {
            res <- NULL
        }
        return(res)
    })

    ## reduce the list to !NULL
    se <- se[!vapply(se, is.null, logical(1))]

    ## Check dimensions
    dims <- sapply(se, dim)
    dims <- c(sum(dims[1,]), max(dims[2,]))
    
    if(any(dim(se) == 0)) {
        index <- GenomicRanges::GRanges()
    }
    else {
        se <- lapply(se, function(y) {
            ## set correct seqinfo
            GenomeInfoDb::seqlevels(rowRanges(y), pruning.mode="coarse") <- GenomeInfoDb::seqlevelsInUse(rowRanges(y))
            return(y)
        })
        slot(settings, "chromosomeList") <- GenomeInfoDb::seqlevels(se[[1]])
        index <- .subsetIndex(se, getIndex(query))
    }

    ## make id slot
    splitid <- .rowRangesFromList(se)
    splitid$id <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(splitid)))
    
    ggd <- new("GenoGAMDataSetList", settings = settings,
               design = design, sizeFactors = sf, index = index,
               ggd = se, id = splitid)
    
    return(ggd)
}

#' @rdname GenoGAMDataSetList-subsetting
setMethod("subsetByOverlaps", c("GenoGAMDataSetList", "GRanges"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"),
                   invert = FALSE, ...) {
              res <- .subsetByOverlapsGDDL(query = query, subject = subject,
                                       maxgap = maxgap, minoverlap = minoverlap,
                                       type = type, invert = invert)
              return(res)
          })

#' @rdname GenoGAMDataSetList-subsetting
setMethod("[", c("GenoGAMDataSetList", "GRanges"), function(x, i) {
    ggd <- subsetByOverlaps(x, i)
    return(ggd)
})

#' @rdname GenoGAMDataSetList-subsetting
setMethod("[[", c("GenoGAMDataSetList", "numeric"), function(x, i) {
    gr <- getIndex(x)[i]
    ggd <- subsetByOverlaps(x,gr)
    return(ggd)
})

## ## Tile computation
## ## ================

## #' Function to retrieve the row coordinates as a list
## #' @param x The GenoGAMDataSet object
## #' @return An integerList with the row numbers for each tile
.getCoordinatesGGDL <- function(x) {

    ## if genome is complete use the fast Bioconductor function
    if(sum(seqlengths(x)) == length(x)) {
        coords <- absoluteRanges(getIndex(x))
    }
    ## otherwise the slower version 'by block'
    else {
        rr <- rowRanges(x)
        coords <- ranges(getIndex(x))
        current <- 1
        for(ii in 1:length(rr)) {
            ov <- IRanges::findOverlaps(rr[[ii]], getIndex(x))
            sh <- S4Vectors::subjectHits(ov)
            qh <- S4Vectors::queryHits(ov)
            l <- range(IRanges::splitAsList(qh, sh))
            l <- IRanges::IRanges(l[,1], l[,2])
            len <- length(l) + current - 1
            end(l) <- end(l) + max(end(coords)[current - 1], 0)
            start(l) <- start(l) + max(end(coords)[current - 1], 0)
            coords[current:len] <- l
            current <- len + 1
        }
    }
    return(coords)
}

## #' Function to establish chunk coordinates
## #' @param x A IRanges object as the output of .getCoordinates
## #' @return The same object as x but with not overlapping ranges
## #' which were cut at the center of the overhang
## .getChunkCoords <- function(x) {
##     if(length(x) == 0) {
##         return(x)
##     }
    
##     start <- c(start(x[1]), ceiling((end(x[-length(x)]) + start(x[-1]))/2))
##     end <- c((start[-1] - 1), end(x[length(x)]))
##     ir <- IRanges(start, end)
##     return(ir)
## }
    
## #' compute metrics for each tile
## #' @param x The GenoGAMDataSet object
## #' @param what A character naming the metric
## #' @param na.rm Should NAs be ignored
## #' @return The metric value
## .MetricsFun <- function(x, what, na.rm = FALSE) {

##     l <- .getCoordinates(x)

##     res <- sapply(colnames(x), function(y) {
##         rle <- IRanges::extractList(assay(x)[[y]], l)
##         eval(call(what, rle, na.rm = na.rm))
##     })
    
##     if(!is(res, "matrix")) {
##         if(is(res, "list") & is.null(dim(res))) {
##             res <- matrix()
##         }
##         else {
##             res <- t(matrix(res))
##         }
##     }

##     colnames(res) <- colnames(x)
##     rownames(res) <- getIndex(x)$id
##     return(res)
## }

## #' Computing metrics
## #'
## #' Computing metrics on each tile of the GenoGAMDataSet object.
## #' All metrics from the Summary generics group, as well as
## #' mean, var, sd, median, mad and IQR are supported.
## #'
## #' @param x A GenoGAMDataSet object
## #' @param ... Additional arguments
## #' @param na.rm Should NAs be dropped. Otherwise the result is NA
## #' @return A matrix with the specified metric computed per tile per column
## #' of the assay data.
## #' @examples
## #' ggd <- makeTestGenoGAMDataSet()
## #' sum(ggd)
## #' min(ggd)
## #' max(ggd)
## #' mean(ggd)
## #' var(ggd)
## #' sd(ggd)
## #' median(ggd)
## #' mad(ggd)
## #' IQR(ggd)
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## #' @rdname GenoGAMDataSet-metrics
## setMethod("Summary", "GenoGAMDataSet", function(x, ..., na.rm = FALSE) {
##     l <- .getCoordinates(x)
        
##     res <- sapply(colnames(x), function(y) {
##         rle <- IRanges::extractList(assay(x)[[y]], l)
##         (getFunction(.Generic))(rle, na.rm = na.rm)
##     })

##     if(!is(res, "matrix")) {
##         if(is(res, "list") & is.null(dim(res))) {
##             res <- matrix()
##         }
##         else {
##             res <- t(matrix(res))
##         }
##     }

##     colnames(res) <- colnames(x)
##     rownames(res) <- getIndex(x)$id
##     return(res)
## })

## #' @rdname GenoGAMDataSet-metrics
## setMethod("mean", "GenoGAMDataSet", function(x) {
##     .MetricsFun(x, "mean")
## })

## #' @rdname GenoGAMDataSet-metrics
## setMethod("var", "GenoGAMDataSet", function(x) {
##     .MetricsFun(x, "var")
## })

## #' @rdname GenoGAMDataSet-metrics
## setMethod("sd", "GenoGAMDataSet", function(x) {
##     .MetricsFun(x, "sd")
## })

## #' @rdname GenoGAMDataSet-metrics
## setMethod("median", "GenoGAMDataSet", function(x) {
##     .MetricsFun(x, "median")
## })

## #' @rdname GenoGAMDataSet-metrics
## setMethod("mad", "GenoGAMDataSet", function(x) {
##     .MetricsFun(x, "mad")
## })

## #' @rdname GenoGAMDataSet-metrics
## setMethod("IQR", "GenoGAMDataSet", function(x) {
##     .MetricsFun(x, "IQR")
## })
