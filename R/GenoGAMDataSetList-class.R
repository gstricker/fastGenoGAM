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
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
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

## .validateChromosomes <- function(object) {
##     cindex <- GenomeInfoDb::seqlevels(object@index)
##     cobject <- GenomeInfoDb::seqlevels(rowRanges(object))
##     if(!all(cindex %in% cobject)) {
##         return("Different chromosomes for data and index objects.")
##     }
##     NULL
## }

## general validate function
.validateGenoGAMDataSetList <- function(object) {
    c(.validateSettingsType(object), .validateDesignType(object),
      .validateSFType(object), .validateIndexType(object),
      ## .validateChromosomes(object),
      .validateGGDType(object),
      .validateIDType(object))
}

S4Vectors::setValidity2("GenoGAMDataSetList", .validateGenoGAMDataSetList)

## Constructor
## ===========

#' GenoGAMDataSetList constructor.
#'
#' GenoGAMDataSetList is the constructor function for the GenoGAMDataSetList-class. 
#' @noRd
GenoGAMDataSetList <- function(...) {
    return(new("GenoGAMDataSetList", ...))
}


#' Make an example /code{GenoGAMDataSet}
#'
#' @param sim Use simulated data (TRUE) or test data from a real experiment
#' @return A /code{GenoGAMDataSet} object
#' @noRd
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

#' @noRd
setMethod("dim", "GenoGAMDataSetList", function(x) {
    dims <- sapply(ggd@ggd, dim)
    res <- c(sum(dims[1,]), max(dims[2,]))
    return(res)
})

#' @noRd
setMethod("length", "GenoGAMDataSetList", function(x) {
    dim(x)[1]
})

#' @noRd
setMethod("seqlengths", "GenoGAMDataSetList", function(x) {
    return(GenomeInfoDb::seqlengths(rowRanges(x)[[1]]))
})

#' @noRd
setMethod("seqlevels", "GenoGAMDataSetList", function(x) {
    return(GenomeInfoDb::seqlevels(rowRanges(x)[[1]]))
})

##' get colData from the first element of the
##' SummarizedExperiment list
##' @noRd
setMethod("colData", "GenoGAMDataSetList", function(x, ...) {
    return(colData(x@ggd[[1]]))
})

##' get colData from the first element of the
##' SummarizedExperiment list
##' @noRd
setMethod("rowRanges", "GenoGAMDataSetList", function(x, ...) {
    lapply(x@ggd, rowRanges)
})


##' get colnames from the first element of the
##' SummarizedExperiment list
##' @noRd
setMethod("colnames", "GenoGAMDataSetList", function(x) {
    return(colnames(x@ggd[[1]]))
})

##' accessor to the index slot
##' @noRd
setMethod("getIndex", "GenoGAMDataSetList", function(object) {
    return(slot(object, "index"))
})

##' The accessor to the list of settings, that were used to generate the tiles.
##' @noRd
setMethod("tileSettings", "GenoGAMDataSetList", function(object) {
    S4Vectors::metadata(getIndex(object))
})

#' make GRanges from list of SummarizedExperiments
#' @noRd
.rowRangesFromList <- function(object) {
    res <- lapply(object, function(x) {
        .extractGR(rowRanges(x))
    })
    return(do.call("c", unname(res)))
}

##' The actual underlying GRanges showing the range of the data.
##' @noRd
setMethod("dataRange", "GenoGAMDataSetList", function(object) {
    res <- .rowRangesFromList(object@ggd)
})

##' A GRanges object representing the chromosomes or chromosome regions
##' on which the model will be computed
##' @noRd
setMethod("getChromosomes", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$chromosomes
})

##' The size of the tiles
##' @noRd
setMethod("getTileSize", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$tileSize
})

##' The size of the chunks
##' @noRd
setMethod("getChunkSize", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$chunkSize
})

##' The size of the overhang (on one side)
##' @noRd
setMethod("getOverhangSize", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$overhangSize
})

##' The total number of tiles
##' @noRd
setMethod("getTileNumber", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$numTiles
})

##' Access to the design slot.
##' @noRd
setMethod("design", "GenoGAMDataSetList", function(object) {
    slot(object, "design")
})

##' Replace method of the design slot.
##' @noRd
setReplaceMethod("design", "GenoGAMDataSetList", function(object, value) {
    newCols <- as.vector(na.omit(.getVars(value)))
    if(!all(newCols %in% colnames(colData(object)))) {
        futile.logger::flog.error("'by' variables could not be found in colData")
        stop("'by' variables could not be found in colData")
    }
    slot(object, "design") <- value
    return(object)
})

##' Access to the sizeFactors slot
##' @noRd
setMethod("sizeFactors", "GenoGAMDataSetList", function(object) {
    sf <- slot(object, "sizeFactors")
    names(sf) <- colnames(object)
    return(sf)
})

##' Replace method of the sizeFactors slot
##' @noRd
setReplaceMethod("sizeFactors", "GenoGAMDataSetList", function(object, value) {
    slot(object, "sizeFactors") <- value
    return(object)
})


## change Settings
## ===============

##' Replace method of the chunkSize parameter,
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

##' Replace method of the tileSize parameter,
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

##' Replace method of the overhangSize parameter,
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

##' Replace method of the tileNumber parameter,
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

## ## Subsetting
## ## ==========
## #' Subset methods for GenoGAMDataSet
## #'
## #' @details
## #' Those are various methods to subset the GenoGAMDataSet object.
## #' By logical statement or GRanges overlap. The '[' subsetter is
## #' just a short version of 'subsetByOverlaps'. The double brackets '[['
## #' offer a subset based on tiles.
## #'
## #' @aliases subsetByOverlaps '[' '[['
## #' @param x,query A GenoGAMDataSet object.
## #' @param subject,i A GRanges object. In case of subsetting by double brackets
## #' 'i' is the index of the tile.
## #' @param maxgap,minoverlap Intervals with a separation of 'maxgap' or
## #' less and a minimum of 'minoverlap' overlapping positions, allowing for
## #' 'maxgap', are considered to be overlapping. 'maxgap' should
## #' be a scalar, non-negative, integer. 'minoverlap' should be a scalar,
## #' positive integer.
## #' @param type By default, any overlap is accepted. By specifying the 'type'
## #' parameter, one can select for specific types of overlap. The types correspond
## #' to operations in Allen's Interval Algebra (see references in). If \code{type}
## #' is \code{start} or \code{end}, the intervals are required to have matching
## #' starts or ends, respectively. While this operation seems trivial, the naive
## #' implementation using \code{outer} would be much less efficient. Specifying
## #' \code{equal} as the type returns the intersection of the \code{start} and
## #' \code{end} matches. If \code{type} is \code{within}, the query interval must
## #' be wholly contained within the subject interval. Note that all matches must
## #' additionally satisfy the \code{minoverlap} constraint described above.
## #'
## #' The \code{maxgap} parameter has special meaning with the special
## #' overlap types. For \code{start}, \code{end}, and \code{equal}, it specifies
## #' the maximum difference in the starts, ends or both, respectively. For
## #' \code{within}, it is the maximum amount by which the query may be wider
## #' than the subject.
## #' @param invert If TRUE, keep only the query ranges that do _not_ overlap
## #' the subject.
## #' @param ... Further arguments. Mostly a logical statement
## #' in case of the 'subset' function. Note that the columnnames
## #' for chromosomes and positions are: 'seqnames' and 'pos'.
## #' @examples
## #' 
## #' # subset by overlaps
## #' ggd <- makeTestGenoGAMDataSet()
## #' SummarizedExperiment::rowRanges(ggd)
## #' gr <- GenomicRanges::GRanges("chrXIV", IRanges(306200,307800))
## #' res <- IRanges::subsetByOverlaps(ggd, gr)
## #' SummarizedExperiment::rowRanges(res)
## #'
## #' # Subset by logical statement
## #' ggd <- makeTestGenoGAMDataSet()
## #' SummarizedExperiment::rowRanges(ggd)
## #' res <- subset(ggd, seqnames == "chrXIV" & pos <= 307000)
## #' SummarizedExperiment::rowRanges(res)
## #' @references
## #' Allen's Interval Algebra: James F. Allen: Maintaining knowledge
## #' about temporal intervals. In: Communications of the ACM.
## #' 26/11/1983. ACM Press. S. 832-843, ISSN 0001-0782
## #' @return A subsetted GenoGAMDataSet object.
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## #' @rdname GenoGAMDataSet-subsetting
## setMethod("subset", "GenoGAMDataSetList", function(x, ...) {
##     if(all(dim(x) == c(0, 0))) return(x)
    
##     settings <- slot(x, "settings")
##     design <- design(x)
##     sf <- sizeFactors(x)

##     ## make ggd slot
##     ## subset all SummarizedExperiments
##     se <- lapply(ggd@ggd, function(se) {
##         res <- subset(se, ...)
##         if(length(res) == 0) {
##             res <- NULL
##         }
##         return(res)
##     })

##     ## reduce the list to !NULL
##     se <- se[!vapply(se, is.null, logical(1))]

##     ## Check dimensions
##     dims <- sapply(se, dim)
##     dims <- c(sum(dims[1,]), max(dims[2,]))

##     ## make index slot
##     ## make empty index if data is empty otherwise trim accordingly
##     if(any(dims == 0)) {
##         index <- GenomicRanges::GRanges()
##     }
##     else {
##         se <- lapply(se, function(y) {
##             ## set correct seqinfo
##             GenomeInfoDb::seqlevels(rowRanges(y), pruning.mode="coarse") <- GenomeInfoDb::seqlevelsInUse(rowRanges(y))
##             return(y)
##         })
##         slot(settings, "chromosomeList") <- GenomeInfoDb::seqlevels(se[[1]])
##         index <- .subsetIndex(se, getIndex(x))
##     }

##     ## make id slot
##     splitid <- .rowRangesFromList(se)
##     splitid$id <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(splitid)))
    
##     ggd <- new("GenoGAMDataSetList", settings = settings,
##                design = design, sizeFactors = sf, index = index,
##                ggd = se, id = splitid)
##     return(ggd)
## })

## #' Method to subset the index of GenoGAMDataSet
## #'
## #' Subsetting the indeces of GenomicTiles based on a SummarizedExperiment.
## #'
## #' @param x A SummarizedExperiment object.
## #' @param index The index of the GenoGAMDataSet object to be subsetted.
## #' @return A GRanges object representing the index
## #' @noRd
## .subsetIndex <- function(se, index) {
##     ## get the data range
##     gpCoords <- .rowRangesFromList(se)
   
##     l <- S4Vectors::metadata(index)
##     l$chromosomes <- gpCoords
##     minWidth <- min(width(gpCoords))
##     if(minWidth < l$tileSize) {
##         l$tileSize <- minWidth
##         l$chunkSize <- minWidth - 2*l$overhangSize
##     }
##     indx <- .makeTiles(l)
##     GenomeInfoDb::seqinfo(indx) <- GenomeInfoDb::seqinfo(gpCoords)
    
##     return(indx)
## }

## .subsetByOverlaps <- function(query, subject, maxgap = 0L, minoverlap = 1L,
##                    type = c("any", "start", "end", "within", "equal"),
##                    invert = FALSE, ...) {
##     if(any((width(subject) %% 2) == 1)) {
##         futile.logger::flog.info("Some subset ranges have odd widths. Rounding to the next even number.")
##         idx <- which((width(subject) %% 2) == 1)
##         width(subject)[idx] <- width(subject)[idx] + 1
##     }          
##     settings <- slot(query, "settings")
##     design <- design(query)
##     sf <- sizeFactors(query)
##     se <- SummarizedExperiment(assays = assays(query),
##                                rowRanges = rowRanges(query),
##                                colData = colData(query))
##     subse <- subsetByOverlaps(se, subject, maxgap = maxgap,
##                               minoverlap = minoverlap,
##                               type=type, invert = invert, ...)
##     if(any(dim(subse) == 0)) {
##         index <- GenomicRanges::GRanges()
##     }
##     else {
##         GenomeInfoDb::seqlevels(rowRanges(subse), pruning.mode="coarse") <- GenomeInfoDb::seqlevelsInUse(rowRanges(subse))
##         slot(settings, "chromosomeList") <- GenomeInfoDb::seqlevels(subse)
##         index <- .subsetIndex(subse, getIndex(query))
##     }
    
##     ggd <- new("GenoGAMDataSet", subse, settings = settings,
##                design = design, sizeFactors = sf, index = index)
##     return(ggd)
## }

## #' @rdname GenoGAMDataSet-subsetting
## setMethod("subsetByOverlaps", c("GenoGAMDataSet", "GRanges"),
##           function(query, subject, maxgap = 0L, minoverlap = 1L,
##                    type = c("any", "start", "end", "within", "equal"),
##                    invert = FALSE, ...) {
##               res <- .subsetByOverlaps(query = query, subject = subject,
##                                        maxgap = maxgap, minoverlap = minoverlap,
##                                        type = type, invert = invert)
##               return(res)
##           })

## #' @rdname GenoGAMDataSet-subsetting
## setMethod("[", c("GenoGAMDataSet", "GRanges"), function(x, i) {
##     ggd <- subsetByOverlaps(x, i)
##     return(ggd)
## })

## #' @rdname GenoGAMDataSet-subsetting
## setMethod("[[", c("GenoGAMDataSet", "numeric"), function(x, i) {
##     gr <- getIndex(x)[i]
##     ggd <- subsetByOverlaps(x,gr)
##     return(ggd)
## })
