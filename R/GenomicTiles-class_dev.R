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
