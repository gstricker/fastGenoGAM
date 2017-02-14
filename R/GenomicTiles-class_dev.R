## #' @rdname getChunkIndex
## #' @export
## setGeneric("getChunkIndex", function(object, ...) standardGeneric("getChunkIndex"))

## #' Compute the index for chunks instead tiles
## #'
## #' The chunk index holds the Granges object that splits the entire dataset in
## #' chunk, that is non-overlapping intervals.
## #'
## #' @rdname getChunkIndex
## #' @param object A /code{GenomicTiles} object.
## #' @param id A vector if tile ids. By default the complete index is returned.
## #' @param ... Additional arguments
## #' @return A /code{GRanges} object representing the index
## #' @examples 
## #' gt <- makeTestGenomicTiles()
## #' getChunkIndex(gt)
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## #' @export
## setMethod("getChunkIndex", "GenomicTiles", function(object, id = NULL) {
##   index <- getIndex(object)
##   if(length(index) == 0) return(index)
##   rowCoords <- getCoordinates(object)
##   gpCoords <- slot(rowRanges(object), "pos_runs")
##   GenomeInfoDb::seqlengths(index) <- rep(NA, length(GenomeInfoDb::seqlengths(index)))
##   index$block <- subjectHits(findOverlaps(index, gpCoords))
    
##   splitIndx <- split(index, index$block)    
  
##   chunkIndex <- lapply(splitIndx, function(y) {
##     start <- c(start(y[1]), ceiling((end(y[-length(y)]) + start(y[-1]))/2))
##     end <- c(ceiling((end(y[-length(y)]) + start(y[-1]))/2 - 1), end(y[length(y)]))
##     start(y) <- start
##     end(y) <- end
##     return(y)
##   })

##   res <- do.call(c, unname(chunkIndex))
##   res$block <- NULL
##   metadata(res) <- metadata(index)
  
##   if(!is.null(id)) res <- res[mcols(res)$id %in% id]
##   return(res)
## })

## ## Index Modification
## ## ===============

## ## The underlying function to unlist indeces, see the genric unlistCoordinates function
## .unlistCoordinates <- function(object, index = NULL, chromosomes = NULL) {
##     if(is.null(index)) index <- getCoordinates(object)
##     if(is.null(chromosomes)) chromosomes <- GenomeInfoDb::seqlevels(index)
    
##     numAssays <- length(assays(object))
##     dims <- dim(object)

##     res <- lapply(chromosomes, function(y) {
##         idx <- index[GenomeInfoDb::seqnames(index) == y]
##         unlistedCoords <- do.call(c, lapply(1:numAssays, function(step) {
##             if(step > 1) shift(ranges(idx), dims[1])
##             else ranges(idx)
##         }))
##         mcols(unlistedCoords)$id <- Rle(factor(rep(idx$id, numAssays)))
##         return(unlistedCoords)
##     })
##     names(res) <- chromosomes
##     return(res)
## }

## setGeneric("unlistCoordinates", function(object, ...) standardGeneric("unlistCoordinates")) 

## setMethod("unlistCoordinates", "GenomicTiles", function(object, index = NULL,
##                                                         chromosomes = NULL) {
##     .unlistCoordinates(object, index = index, chromosomes = chromosomes)
## })

## setGeneric("unlistIndexCoordinates",
##            function(object, ...) standardGeneric("unlistIndexCoordinates")) 

## setMethod("unlistIndexCoordinates", "GenomicTiles", function(object, chromosomes = NULL) {
##     rindex <- getIndexCoordinates(object)
##     .unlistCoordinates(object, chromosomes = chromosomes, index = rindex)
## })

## setGeneric("unlistIndex", function(object, ...) standardGeneric("unlistIndex")) 

## setMethod("unlistIndex", "GenomicTiles", function(object, chromosomes = NULL) {
##     tindex <- getIndex(object)
##     .unlistCoordinates(object, chromosomes = chromosomes, index = tindex)
## })







## #' @rdname GenomicTiles-view
## #' @export
## setGeneric("view", function(object, ...) standardGeneric("view"))

## #' View the dataset
## #'
## #' Cbinding the columns all together and coercing to data.frame
## #'
## #' @param object A \code{GenomicTiles} object
## #' @param ranges A \code{GRanges} object. Makes it possible to
## #' select regions by \code{GRanges}. Either ranges or seqnames, start and
## #' end must be supplied
## #' @param seqnames A chromosomes name. Either ranges or seqnames, start and
## #' end must be supplied
## #' @param start A start site. Either ranges or seqnames, start and
## #' end must be supplied
## #' @param end An end site. Either ranges or seqnames, start and
## #' end must be supplied
## #' @param ... Additional arguments
## #' @return A data.frame of the selected data.
## #' @examples
## #' gt <- makeTestGenomicTiles()
## #' gr <- GRanges(c("chrI", "chrII"), IRanges(c(1, 10), c(40, 30)))
## #' head(view(gt, ranges = gr))
## #' head(view(gt, seqnames = "chrI", start = 1, end = 20))
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## #' @rdname GenomicTiles-view
## #' @export
## setMethod("view", c("GenomicTiles"), function(object, ranges = NULL, seqnames = NULL,
##                                               start = NULL, end = NULL) {
##     if((is.null(seqnames) | is.null(start) | is.null(end)) & is.null(ranges)) {
##         return(as.data.frame(object))
##     }
##     if(!is.null(ranges)) {
##         res <- .exactSubsetByOverlaps(object, ranges)
##     }
##     else {
##         res <- subset(object, seqnames == seqnames & pos >= start & pos <= end)
##     }
##     return(as.data.frame(res))
## })

## ## Tile computation
## ## ================

## #' compute metrics for each tile
## .MetricsFun <- function(x, what, na.rm = FALSE) {

##     indx <- getIndexCoordinates(x)
    
##     res <- lapply(1:length(assays(x)), function(ii) {
##         sapply(rownames(colData(x)), function(y) {
##             rle <- extractList(assay(x, ii)[[y]], ranges(indx))
##             eval(call(what, rle, na.rm = na.rm))
##         })
##     })

##     names(res) <- names(assays(x))
##     return(res)
## }

## #' Computing metrics
## #'
## #' Computing metrics on each tile of the \code{GenomicTiles} object.
## #' So far all metrics from the Summary generics group, as well as
## #' mean, var, sd, median, mad and IQR are supported.
## #'
## #' @param x A \code{GenomicTiles} object
## #' @param ... Additional arguments
## #' @param na.rm Should NAs be dropped. Otherwise the result is NA
## #' @return A list of as many elements as there are assays.
## #' Each element contains of a matrix with the specified
## #' metric computed per tile per column of the assay data.
## #' @examples
## #' gt <- makeTestGenomicTiles()
## #' sum(gt)
## #' min(gt)
## #' max(gt)
## #' mean(gt)
## #' var(gt)
## #' sd(gt)
## #' median(gt)
## #' mad(gt)
## #' IQR(gt)
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## #' @rdname GenomicTiles-metrics
## setMethod("Summary", "GenomicTiles", function(x, ..., na.rm = FALSE) {
##     indx <- getIndexCoordinates(x)
    
##     res <- lapply(1:length(assays(x)), function(ii) {
##         sapply(rownames(colData(x)), function(y) {
##             rle <- extractList(assay(x, ii)[[y]], ranges(indx))
##             (getFunction(.Generic))(rle, na.rm = na.rm)
##         })
##     })

##     names(res) <- names(assays(x))
##     return(res)
## })

## #' @rdname GenomicTiles-metrics
## setMethod("mean", "GenomicTiles", function(x) {
##     .MetricsFun(x, "mean")
## })

## #' @rdname GenomicTiles-metrics
## setMethod("var", "GenomicTiles", function(x) {
##     .MetricsFun(x, "var")
## })

## #' @rdname GenomicTiles-metrics
## setMethod("sd", "GenomicTiles", function(x) {
##     .MetricsFun(x, "sd")
## })

## #' @rdname GenomicTiles-metrics
## setMethod("median", "GenomicTiles", function(x) {
##     .MetricsFun(x, "median")
## })

## #' @rdname GenomicTiles-metrics
## setMethod("mad", "GenomicTiles", function(x) {
##     .MetricsFun(x, "mad")
## })

## #' @rdname GenomicTiles-metrics
## setMethod("IQR", "GenomicTiles", function(x) {
##     .MetricsFun(x, "IQR")
## })

## ## Cosmetics
## ## =========

## .showGenomicTiles <- function(object) {
##     cl <- class(object)
##     dims <- dim(object)
##     as <- assays(object)
##     rd <- names(rowData(object))
##     md <- unique(names(colData(object)))
##     cnames <- colnames(object)
##     cdata <- names(colData(object))

##     if(length(tileSettings(object)) != 0) {
##         tsize <- tileSettings(object)$tileSize
##         tname <- "tiles"
##         chunk <- tileSettings(object)$chunk
##         if (chunk) {
##             tsize <- tileSettings(object)$chunkSize
##             tname <- "chunks"
##         }
##         unit <- ""
##         if(!is.null(tsize)) {
##             unit = "bp"
##             if(tsize/1000 > 1) {
##                 unit <- "kbp"
##                 tsize <- tsize/1000
##             }
##             if(tsize/1e6 > 1) {
##                 unit <- "Mbp"
##                 tsize <- tsize/1e6
##             }
##         }
##         chroms <- GenomeInfoDb::seqlevels(tileSettings(object)$chromosomes)
##     }
##     tnum <- length(getIndex(object))
    
    
##     cat("class:", cl, "\n")
##     cat("dimension:", dims, "\n")
##     cat(paste0("assays(", length(as), "):"), names(as), "\n")
##     cat(paste0("position variables(", length(rd), "):"), rd, "\n")
##     cat(paste0("sample variables(", length(md), "):"), md, "\n")
##     cat(paste0("samples(", length(cnames), "):"), cnames, "\n")
##     if(length(tileSettings(object)) != 0) {
##         cat(paste0(tname, " size: ", tsize, unit), "\n")
##         cat(paste0("number of ", tname, ": ", tnum), "\n")
##         cat("chromosomes:", chroms, "\n")
##     }
## }

## ## Show method for GenomicTiles.
## setMethod("show", "GenomicTiles", function(object) {
##     .showGenomicTiles(object)
## })


