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


## ## Subsetting
## ## ==========

## #' Method to subset the indeces of GenomicTiles.
## #'
## #' Subsetting the indeces of GenomicTiles based on a SummarizedExperiment.
## #'
## #' @param x A SummarizedExperiment object.
## #' @param index The index of a GenomicTiles to be subsetted.
## #' @return A list of two GRanges objects: the subsetted index and the
## #' subsetted coordinates.
## .subsetIndeces <- function(se, index) {
##     gpCoords <- rowRanges(se)@pos_runs
##     gpDims <- length(gpCoords)
##     ir <- IRanges(start = cumsum(c(1, end(gpCoords)[-gpDims[1]])),
##                   end = cumsum(width(gpCoords)))
##     coords <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(gpCoords), ir)
##     l <- metadata(index)
##     l$chromosomes <- gpCoords
    
##     indx <- .makeTiles(l)
##     coords <- .makeCoordinates(GPos(coords))
##     GenomeInfoDb::seqinfo(indx) <- GenomeInfoDb::seqinfo(gpCoords)
##     GenomeInfoDb::seqlevels(coords) <- GenomeInfoDb::seqlevels(coords)[GenomeInfoDb::seqlevels(coords) %in%
##                                            unique(GenomeInfoDb::seqnames(coords))]
##     GenomeInfoDb::seqlevels(indx) <- GenomeInfoDb::seqlevels(indx)[GenomeInfoDb::seqlevels(indx) %in%
##                                            unique(GenomeInfoDb::seqnames(indx))]
##     return(list(index = indx, coordinates = coords))
## }

## #' Method to subset GenomicTiles.
## #'
## #' Subsetting the GenomicTiles based on the subset method of
## #' SummarizedExperiment. Additionally the index and coordinates gets
## #' subsetted.
## #'
## #' @param x A GenomicTiles object.
## #' @param ... Any subset option accepted by SummarizedExperiment. Usually
## #' in accordance to base::subset.
## #' @return A subsetted GenomicTiles object.
## .subsetGenomicTiles <- function(x,...){
##     if(all(dim(x) == c(0, 0))) return(x)
##     se <- subset(SummarizedExperiment(assays = assays(x),
##                                       rowRanges = rowRanges(x),
##                                       colData = colData(x)), ...)
##     GenomeInfoDb::seqlevels(rowRanges(se), force = TRUE) <- GenomeInfoDb::seqlevelsInUse(rowRanges(se))

##     indeces <- .subsetIndeces(se, getIndex(x))
##     GenomeInfoDb::seqlevels(metadata(indeces$index)$chromosomes, force = TRUE) <- GenomeInfoDb::seqlevelsInUse(rowRanges(se))
##     return(new("GenomicTiles", index = indeces$index,
##                coordinates = indeces$coordinates, se))
## }


## #' Subset method for /code{GenomciTiles}
## #'
## #' Subsetting the /code{GenomicTiles} by a logical statement
## #'
## #' @param x A /code{GenomicTiles} object.
## #' @param ... Further arguments. Mostly a logical statement.
## #' Note that the columnnames for chromosomes and positions
## #' are: seqnames and pos.
## #' @return A subsetted /code{GenomicTiles} object.
## #' @examples
## #' gt <- makeTestGenomicTiles()
## #' res <- subset(gt, seqnames == "chrI" & pos <= 50)
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## setMethod("subset", "GenomicTiles", function(x, ...) {
##     .subsetGenomicTiles(x, ...)
## })

## #' Method to subset GenomicTiles by a GRanges object.
## #'
## #' Subsetting the GenomicTiles by a GRanges object based on the
## #' subsetByOverlaps method applying to SummarizedExperiment.
## #' Additionally the index and coordinates gets subsetted.
## #'
## #' @param query A GenomicTiles object.
## #' @param subject A GRanges object.
## #' @param ... Any other parameters applicable to subset a
## #' SummarizedExperiment by overlaps.
## #' @return A subsetted GenomicTiles object.
## .subsetGenomicTilesByOverlaps <- function(query, subject, ...){
##     se <- subsetByOverlaps(SummarizedExperiment(assays = assays(query),
##                                       rowRanges = rowRanges(query)),
##                            subject, ...)
##     indeces <- .subsetIndeces(se, getIndex(query))
##     return(new("GenomicTiles", index = indeces$index,
##                coordinates = indeces$coordinates, se))
## }

## .subsetByOverlaps <- function(query, subject, maxgap, minoverlap,
##                               type, ...) {
##     rowRanges <- rowRanges(query)
##     assay <- assays(query)
##     index <- getIndex(query)
##     colData <- colData(query)

##     tiles <- subsetByOverlaps(index, subject)
##     se <- subsetByOverlaps(SummarizedExperiment(assay, rowRanges = rowRanges, colData = colData),
##                            tiles, ...)
##     GenomeInfoDb::seqlevels(rowRanges(se), force = TRUE) <- GenomeInfoDb::seqlevelsInUse(rowRanges(se))
##     GenomeInfoDb::seqlevels(tiles, force = TRUE) <- GenomeInfoDb::seqlevelsInUse(tiles)
##     if(length(metadata(tiles)) > 0) {
##         metadata(tiles)$numTiles <- length(tiles)
##         GenomeInfoDb::seqlevels(metadata(tiles)$chromosomes, force = TRUE) <- GenomeInfoDb::seqlevelsInUse(tiles)
##     }
    
##     coords <- .makeCoordinates(rowRanges(se))
    
##     res <- new("GenomicTiles", index = tiles,
##                coordinates = coords, se)
##     return(res)
## }

## .exactSubsetByOverlaps <- function(query, subject, ...) {
##     rowRanges <- rowRanges(query)
##     assay <- assays(query)
##     settings <- tileSettings(query)
##     colData <- colData(query)

##     se <- subsetByOverlaps(SummarizedExperiment(assay, rowRanges = rowRanges, colData = colData),
##                            subject, ...)
##     GenomeInfoDb::seqlevels(rowRanges(se), force = TRUE) <- GenomeInfoDb::seqlevelsInUse(rowRanges(se))
    
##     settings$chromosomes <- slot(rowRanges(se), "pos_runs")
##     tiles <- .makeTiles(settings)
##     GenomeInfoDb::seqlevels(tiles, force = TRUE) <- GenomeInfoDb::seqlevelsInUse(rowRanges(se))
##     if(length(metadata(tiles)) > 0) {
##         GenomeInfoDb::seqlevels(metadata(tiles)$chromosomes, force = TRUE) <- GenomeInfoDb::seqlevelsInUse(tiles)
##     }
    
##     coords <- .makeCoordinates(rowRanges(se))
    
##     res <- new("GenomicTiles", index = tiles,
##                coordinates = coords, se)
##     return(res)
## }

## #' Subset by overlaps method for \code{GenomciTiles}
## #'
## #' Subsetting the \code{GenomicTiles} by a \code{GRanges} object
## #'
## #' @param query A \code{GenomicTiles} object.
## #' @param subject A \code{GRanges} object
## #' @param maxgap,minoverlap Intervals with a separation of \code{maxgap} or
## #' less and a minimum of \code{minoverlap} overlapping positions, allowing for
## #' \code{maxgap}, are considered to be overlapping.  \code{maxgap} should
## #' be a scalar, non-negative, integer. \code{minoverlap} should be a scalar,
## #' positive integer.
## #' @param type By default, any overlap is accepted. By specifying the \code{type}
## #' parameter, one can select for specific types of overlap. The types correspond
## #' to operations in Allen's Interval Algebra (see references). If \code{type}
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
## #' @param ... Additional parameters
## #' @return A subsetted \code{GenomicTiles} object.
## #' @examples
## #' gt <- makeTestGenomicTiles()
## #' gr <- GRanges(c("chrI", "chrII"), IRanges(c(1, 120), c(40, 150)))
## #' res <- subsetByOverlaps(gt, gr)
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## setMethod("subsetByOverlaps", c("GenomicTiles", "GRanges"),
##           function(query, subject, maxgap=0L, minoverlap=1L,
##                       type=c("any", "start", "end", "within", "equal"), ...) {
##               .subsetByOverlaps(query, subject, maxgap = maxgap,
##                                 minoverlap = minoverlap,
##                                 type = type, ...)
##           })

## #' Method to extract all GenomicTiles at once.
## #'
## #' Extracting and converting all GenomicTiles at once using the complete
## #' assays. This is more memory consuming, but faster for smaller datasets.
## #'
## #' @param gt A GenomicTiles object.
## #' @param index A GRanges object.
## #' @param chromosomes A character vector of chromosome names.
## #' @return A SimpleList of DataFrames with as many elements
## #' as there are tiles.
## .extractFullGenomicTiles <- function(gt, index, chromosomes) {
##     rowindx <- getIndexCoordinates(gt, index = index)
##     coords <- unlistCoordinates(gt, index = rowindx,
##                                 chromosomes = chromosomes)
##     fullIndex <- sort(do.call(c, unname(coords)))
##     df <- as(gt, "DataFrame")
##     meta <- metadata(df)
##     vec <- rep(mcols(fullIndex)$id, width(fullIndex))

##     applyVec <- IRanges(start=(0:(length(fullIndex)/length(rowindx) - 1))*length(rowindx) + 1,
##                         end = (1:(length(fullIndex)/length(rowindx)))*length(rowindx))

##     ## From in to out: For each assay make DataFrameList,
##     ## combine and unlist the List to get a DataFrame of tiles that are consecutively placed.
##     ## Because of overlaps this DataFrame is bigger than the original and the same size
##     ## as "vec".
##     dfList <- splitAsList(unlist(do.call(c, lapply(1:length(applyVec), function(y) {
##         extractList(df, fullIndex[start(applyVec)[y] : end(applyVec)[y]])
##     }))), vec)

##     metadata(dfList) <- meta
##     return(dfList)
## }

## #' Method to extract all GenomicTiles by blocks.
## #'
## #' Extracting and converting all GenomicTiles iteratively by blocks of
## #' GenomicTiles. This is less memory consuming, but slower for smaller
## #' datasets.
## #'
## #' @param gt A GenomicTiles object.
## #' @param index A GRanges object.
## #' @param blocks An integer value representing the number of blocks.
## #' @return A SimpleList of DataFrames with as many elements
## #' as there are tiles.
## .extractGenomicTilesByBlocks <- function(gt, index, blocks) {
##     blockSize <- length(index)/blocks
##     blockList <- unname(split(mcols(index)$id,
##                       ceiling(seq_along(mcols(index)$id)/blockSize)))
    
##     dfList <- do.call(c, bplapply(blockList, function(y) {
##         suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))
##         blockindx <- index[match(y, index$id),]
##         subgt <- .subsetByOverlaps(gt, blockindx)
##         df <- DataFrame(subgt)
##         coords <- unlistCoordinates(subgt,
##                                     getIndexCoordinates(subgt))
##         fullIndex <- do.call(c, unname(coords))
##         vec <- rep(mcols(fullIndex)$id, width(fullIndex))
##         dfList <- splitAsList(df, vec)
##         return(dfList)
##     }))
##     return(dfList)
## }

## .extractGenomicTilesByIndex <- function(gt, index, size = 3e9) {
##     chromosomes <- GenomeInfoDb::seqlevels(index)
##     gtDims <- dim(gt)
##     if(all(gtDims == c(0, 0))) return(DataFrameList())
##     space <- getChunkSize(gt)*length(index)
    
##     if(space <= size) {
##         dfList <- .extractFullGenomicTiles(gt, index, chromosomes)
##     }
##     else {
##         dfList <- .extractGenomicTilesByBlocks(gt, index,
##                                                ceiling(space/size))
##     }
##     return(dfList)
## }

## #' Get a tile.
## #'
## #' Extracting one tile from a GenomicTiles object.
## #'
## #' @param object A GenomicTiles Object.
## #' @param id A tile id as an integer
## #' @return A data.frame of the tile.
## .getTile <- function(gt, id) {
##     index <- getIndex(gt)
##     range <- ranges(index[index$id == id,])
##     chrom <- as.character(GenomeInfoDb::seqnames(index[index$id == id,]))
##     gtSubset <- subset(gt, GenomeInfoDb::seqnames(gt) == chrom &
##                        pos >= start(range) & pos <= end(range))

##     df <- DataFrame(gtSubset)
##     ##names(df)[1] <- "chromosome"
##     ## df <- melt(df, measure.vars = rownames(colData(gtSubset)), variable.name = "sample")
##     return(DataFrameList(df))
## }

## #' @rdname getTile
## #' @export
## setGeneric("getTile", function(object, id, ...) standardGeneric("getTile"))

## #' Tile extraction as a DataFrame
## #'
## #' Extracting one or multiple tiles from a \code{GenomicTiles} object
## #' and coercing them to a DataFrameList.
## #'
## #' @param object A \code{GenomicTiles} object
## #' @param id A vector of tile ids
## #' @param size The maximal number of rows that should be handled at once.
## #' If the dataset is bigger it will be processed in chunks. This is to lower
## #' memory consumption on big datasets, which in turn is slower.
## #' @param ... Additional arguments
## #' @return A \code{SimpleDataFrameList}
## #' @examples
## #' gt <- makeTestGenomicTiles()
## #' getTile(gt, 1:3)
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## #' @rdname getTile
## #' @export
## setMethod("getTile", "GenomicTiles", function(object, id, size = 3e9) {
##     if(missing(id)) {
##         indx <- getIndex(object)
##         id <- mcols(indx)$id
##         res <- .extractGenomicTilesByIndex(object,
##                                            indx,
##                                            size = size)
##     }
##     else {
##         if(length(id) == 1) res <- .getTile(object, id)
##         if(length(id) > 1) {
##             indx <- getIndex(object, id)
##             res <- .extractGenomicTilesByIndex(object, indx, size = size)
##         }
##     }
##     meta <- tileSettings(object)
##     meta$chunks <- getChunkIndex(object, id)
##     metadata(res) <- meta
##     return(res)
## })

## .subsetByDoubleBrackets <- function(x, i, j) {
##     if(!missing(j)) {
##         stop("Wrong number of dimensions")
##     }

##     if(length(i) > 1L) {
##         stop("attempt to extract more than one element")
##     }

##     return(.getTile(x, i))
## }

## #' Providing pseudo-list functionality
## #'
## #' Getting a specific tile
## #'
## #' @param x A \code{GenomicTiles} object
## #' @param i An integer (for '[[') or a \code{GRanges} object (for '[')
## #' @return A \code{DataFrame} (for '[[') or a subsetted \code{GenomicTiles} object (for '[')
## #' @rdname GenomicTiles-brackets
## setMethod("[[", c("GenomicTiles", "numeric"),
##           function(x, i) .subsetByDoubleBrackets(x, i))

## #' @rdname GenomicTiles-brackets
## setMethod("[", c("GenomicTiles", "GRanges"),
##           function(x, i) .exactSubsetByOverlaps(x, i))

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


