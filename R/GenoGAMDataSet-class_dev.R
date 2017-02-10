## ====================
## GenoGAMDataset class
## ====================
#' @include GenoGAMSettings-class.R
NULL

#' GenoGAMDataSet
#'
#' The GenoGAMDataSet class contains the pre-processed raw data and
#' additional slots that define the input and framework for the model.
#' It extends the RangedSummarizedExperiment class by adding an index
#' that defines ranges on the entire genome, mostly for purposes of
#' parallel evaluation. Furthermore adding a couple more slots to hold
#' information such as experiment design. It also contains the
#' \code{\link[fastGenoGAM]{GenoGAMSettings}} class that defines global
#' settings for the session. For information on the slots inherited from
#' SummarizedExperiment check the respective class.
#' 
#' @slot settings The global and local settings that will be used to compute the
#' model.
#' @slot design The formula describing how to evaluate the data. See details.
#' @slot sizeFactors The normalized values for each sample. A named numeric vector.
#' @slot index A GRanges object representing an index of the ranges defined 
#' on the genome. Mostly used to store tiles.
#' @name GenoGAMDataSet-class
#' @rdname GenoGAMDataSet-class
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @exportClass GenoGAMDataSet
setClass("GenoGAMDataSet",
         contains = "RangedSummarizedExperiment",
         slots = list(settings = "GenoGAMSettings",
                      design = "formula", sizeFactors = "numeric",
                      index = "GRanges"),
         prototype = list(settings = GenoGAMSettings(),
                          design = ~ 1, sizeFactors = numeric(), 
                          index = GenomicRanges::GRanges()))

## Validity
## ========

.validateSettingsType <- function(object) {
    if(class(slot(object, "settings")) != "GenoGAMSettings") {
        return("'settings' must be a GenoGAMSettings object")
    }
    NULL
}

.validateDesignType <- function(object) {
    if(class(slot(object, "design")) != "formula") {
        return("'design' must be a formula object")
    }
    NULL
}

.validateSFType <- function(object) {
    if(class(slot(object, "sizeFactors")) != "numeric") {
        return("'sizeFactors' must be a numeric object")
    }
    NULL
}

#' Validating the correct type
.validateIndexType <- function(object) {
    if(class(object@index) != "GRanges") {
        return("'index' must be a GRanges object")
    }
    NULL
}

.validateChromosomes <- function(object) {
    cindex <- GenomeInfoDb::seqlevels(object@index)
    cobject <- GenomeInfoDb::seqlevels(rowRanges(object))
    if(!all(cindex %in% cobject)) {
        return("Different chromosomes for data and index objects.")
    }
    NULL
}

## general validate function
.validateGenoGAMDataSet <- function(object) {
    c(.validateSettingsType(object), .validateDesignType(object),
      .validateSFType(object), .validateIndexType(object),
      .validateChromosomes(object))
}

setValidity2("GenoGAMDataSet", .validateGenoGAMDataSet)


## Constructor
## ===========

#' GenoGAMDataSet constructor.
#'
#' GenoGAMDataSet is the constructor function for the GenoGAMDataSet-class. 
#'
#' @aliases getIndex tileSettings dataRange getChromosomes getTileSize getChunkSize getOverhangSize getTileNumber design sizeFactors
#' @param experimentDesign Either a character object specifying the path to a
#' delimited text file (the delimiter will be determined automatically),
#' a data.frame specifying the experiment design or a RangedSummarizedExperiment
#' object with the GPos class being the rowRanges. See details for the structure
#' of the experimentDesign.
#' @param chunkSize An integer specifying the size of one chunk in bp.
#' @param overhangSize An integer specifying the size of the overhang in bp.
#' As the overhang is taken to be symmetrical, only the overhang of one side
#' should be provided.
#' @param design A formula object. See details for its structure.
#' @param directory The directory from which to read the data. By default
#' the current working directory is taken.
#' @param settings A GenoGAMSettings object. Not needed by default, but might
#' be of use if only specific regions should be read in
#' See \code{\link{GenoGAMSettings}}.
#' @param hdf5 Should the data be stored on HDD in HDF5 format? By default this
#' is disabled, as the Rle representation of count data already provides a
#' decent compression of the data. However in case of large organisms, a complex
#' experiment design or just limited memory, this might further decrease the
#' memory footprint. Note this only applies to the input count data, results are
#' usually stored in HDF5 format due to their space requirements for type double.
#' @param ... Further parameters, mostly for arguments of custom processing
#' functions or to specify a different method for fragment size estimation.
#' See details for further information.
#' @param object For use of S4 methods. The GenoGAMDataSet object.
#' @param value For use of S4 methods. The value to be assigned to the slot.
#' @return An object of class \code{\link{GenoGAMDataSet}} or the respective slot.
#' @section Config:
#' 
#' The config file/data.frame contains the actual experiment design. It must
#' contain at least three columns with fixed names: 'ID', 'file' and 'paired'.
#'
#' The field 'ID' stores a unique identifier for each alignment file.
#' It is recommended to use short and easy to understand identifiers because
#' they are subsequently used for labelling data and plots.
#'
#' The field 'file' stores the BAM file name.
#'
#' The field 'paired', values TRUE for paired-end sequencing data, and FALSE for
#' single-end sequencing data.
#'
#' All other columns are stored in the colData slot of the GenoGAMDataSet
#' object. Note that all columns which will be used for analysis must have at
#' most two conditions, which are for now restricted to 0 and 1. For example,
#' if the IP data schould be corrected for input, then the input will be 0
#' and IP will be 1, since we are interested in the corrected IP. See examples.
#'
#' @section Design/Formula:
#' 
#' Design must be a formula. At the moment only the following is
#' possible: Either ~ s(x) for a smooth fit over the entire data or
#' s(x, by = "myColumn"), where 'myColumn' is a column name
#' in the experimentDesign. Any combination of this is possible:
#' 
#' ~ s(x) + s(x, by = "myColumn") + s(x, by = ...) + ...
#' 
#' For example the formula for correcting IP for input would look like this:
#' 
#' ~ s(x) + s(x, by = "experiment")
#'
#' where 'experiment' is a column with 0s and 1s, with the ip samples annotated
#' with 1 and input samples with 0.
#''
#' @section Further parameters:
#' 
#' In case of single-end data it might be usefull to specify a different
#' method for fragment size estimation. The argument 'shiftMethod' can be
#' supplied with the values 'coverage' (default), 'correlation' or 'SISSR'.
#' See ?chipseq::estimate.mean.fraglen for explanation.
#' @examples
#' \dontrun{
#' myConfig <- data.frame(ID = c("input","ip"),
#'                   file = c("myInput.bam", "myIP.bam"),
#'                   paired = c(FALSE, FALSE),
#'                   experiment = factor(c(0,1)),
#'                   stringsAsFactors = FALSE) 
#' myConfig2 <- data.frame(ID = c("wildtype1","wildtype2",
#'                               "mutant1", "mutant2"),
#'                   file = c("myWT1.bam", "myWT2.bam"
#'                            "myMutant1.bam", "myMutant2.bam"),
#'                   paired = c(FALSE, FALSE, FALSE, FALSE),
#'                   experiment = factor(c(0, 0, 1, 1)),
#'                   stringsAsFactors = FALSE)
#' 
#' gdd <- GenoGAMDataSet(myConfig, chunkSize = 2000,
#' overhangSize = 250, design = ~ s(x) + s(x, by = "experiment")
#' ggd <- GenoGAMDataSet(myConfig2, chunkSize = 2000,
#' overhangSize = 250, design = ~ s(x) + s(x, by = "experiment"))
#' }
#'
#' ## build from SummarizedExperiment
#' gr <- GPos(GRanges("chr1", IRanges(1, 10000)))
#' seqlengths(gr) <- 1e6
#' df <- DataFrame(colA = 1:10000, colB = round(runif(10000)))
#' se <- SummarizedExperiment(rowRanges = gr, assays = list(df))
#' ggd <- GenoGAMDataSet(se, chunkSize = 2000, overhangSize = 250, 
#'                       design = ~ s(x) + s(x, by = "experiment"))
#' ggd
#' @name GenoGAMDataSet
#' @rdname GenoGAMDataSet-class
#' @export
GenoGAMDataSet <- function(experimentDesign, chunkSize, overhangSize, design,
                           directory = ".", settings = NULL, hdf5 = FALSE, ...) {

    futile.logger::flog.info("Creating GenoGAMDataSet")

    if(missing(experimentDesign)) {
        futile.logger::flog.debug("No input provided. Creating empty GenoGAMDataSet")
        return(new("GenoGAMDataSet"))
    }

    input <- paste0("Building GenoGAMDataSet with the following parameters:\n",
                    "  Class of experimentDesign: ", class(experimentDesign), "\n",
                    "  Chunk size: ", chunkSize, "\n",
                    "  Overhang size: ", overhangSize, "\n",
                    "  Design: ", design, "\n",
                    "  Directory: ", directory, "\n",
                    "  Settings: ", settings, "\n",
                    "  HDF5: ", hdf5, "\n",
                    "  Further parameters: ", list(...), "\n")
    futile.logger::flog.debug(input)

    if(is.null(settings)) {
        settings <- GenoGAMSettings()
    }

    if(class(experimentDesign) == "RangedSummarizedExperiment") {
        futile.logger::flog.trace("Building GenoGAMDataSet from SummarizedExperiment object")
        gt <- .GenoGAMDataSetFromSE(se = experimentDesign,
                                    chunkSize = chunkSize,
                                    overhangSize = overhangSize,
                                    design = design,
                                    settings = settings, ...)
    }
    else {
        futile.logger::flog.trace(paste0("Building GenoGAMDataSet from config file: ", experimentDesign))
        gt <- .GenoGAMDataSetFromConfig(config = experimentDesign,
                                        chunkSize = chunkSize,
                                        overhangSize = overhangSize,
                                        design = design,
                                        directory = directory,
                                        settings = settings, 
                                        hdf5 = hdf5, ...)
    }
    
    return(gt)
}

#' The underlying function to build a GenoGAMDataSet from a
#' SummarizedExperiment
#' 
#' @noRd
.GenoGAMDataSetFromSE <- function(se, chunkSize, overhangSize,
                                    design, settings, ...) {

    ## make tiles
    gr <- rowRanges(se)@pos_runs

    ## check for overlapping ranges
    if(sum(countOverlaps(gr)) > length(gr)) {
        stop("Overlapping regions encountered. Please reduce ranges and data first.")
    }

    if(any(is.na(seqlengths(se)))) {
        stop("Sequence lengths missing in the Seqinfo object of SummarizedExperiment")
    }

    l <- list(chromosomes = gr,
              chunkSize = chunkSize,
              overhangSize = overhangSize)
    tiles <- .makeTiles(l)

    ## initiate size factors
    sf <- rep(0, length(colnames(se)))
    names(sf) <- colnames(se)

    ## update chromosome list
    if(is.null(slot(settings, "chromosomeList"))) {
        slot(settings, "chromosomeList") <- GenomeInfoDb::seqlevels(gr)
    }

    ggd <- new("GenoGAMDataSet", se, settings = settings,
               design = design, sizeFactors = sf, index = tiles)


    ## check if everything was set fine
    correct <- checkObject(ggd)
    if(!correct) stop()
    
    return(ggd)   
}

#' A function to produce a GRanges index from a list of settings.
#'
#' @param l A list of settings involving:
#' chunkSize, chromosomes, tileSize and overhangSize
#' @return A GRanges object of the tiles.
#'
#' @noRd
.makeTiles <- function(l) {

    if(length(l) == 0) return(GenomicRanges::GRanges())
    if(l$overhangSize < 0) stop("Overhang size must be equal or greater than 0")
    if(l$chunkSize < 1000) stop("Chunk size must be equal or greater than 1000")
    if(length(l$chromosomes) == 0) stop("Chromosome list should contain at least one entry")

    input <- paste0("Building Tiles with the following parameters:\n",
                    "  Chunk size: ", l$chunkSize, "\n",
                    "  Overhang size: ", l$overhangSize, "\n",
                    "  Chromosomes: ", l$chromosomes, "\n")
    futile.logger::flog.debug(input)

    l$tileSize <- l$chunkSize + 2*l$overhangSize
    futile.logger::flog.trace(paste0("GenoGAMDataSet: Tile size computed to be ", l$tileSize))

    ## deal with overlapping ranges to reduce complexity and redundancy
    l$chromosomes <- reduce(l$chromosomes)
    lambdaFun <- function(y, sl) {

        ## change to GenoGAM as soon as ready
        suppressPackageStartupMessages(require(fastGenoGAM, quietly = TRUE))

        ## generate break points for chunks
        nchunks <- ceiling(IRanges::width(y)/sl$chunkSize)
      
        startSeq <- seq(0, nchunks - 1, length.out = nchunks) * sl$chunkSize
        endSeq <- seq(1, nchunks - 1, length.out = (nchunks - 1)) * sl$chunkSize - 1
        starts <- IRanges::start(y) + startSeq
        ends <- c(IRanges::start(y) + endSeq, IRanges::end(y))
        ir <- IRanges::IRanges(starts, ends)
        chunks <- GenomicRanges::GRanges(seqnames = seqnames(y), ir)
        seqinfo(chunks) <- seqinfo(y)

        ## flank to tiles
        if(sl$tileSize == sl$chunkSize) {
            tiles <- chunks
        }
        else {
            ov <- round(IRanges::width(chunks)/2)
            centeredChunks <- suppressWarnings(IRanges::shift(chunks, ov))
            ir <- suppressWarnings(IRanges::flank(centeredChunks, round(sl$tileSize/2), both = TRUE))
            tiles <- suppressWarnings(IRanges::trim(ir))
        }

        ## adjust first tile
        startsToResize <- which(IRanges::start(tiles) < IRanges::start(y))
        trimmedEnd <- IRanges::end(tiles[startsToResize]) + IRanges::start(y) - IRanges::start(tiles[startsToResize])
        IRanges::end(tiles[startsToResize]) <- min(trimmedEnd, IRanges::end(y))
        IRanges::start(tiles[startsToResize]) <- IRanges::start(y)

        ## adjust last tile
        endsToResize <- which(IRanges::end(tiles) > IRanges::end(y))
        trimmedStarts <- IRanges::start(tiles[endsToResize]) - IRanges::end(tiles[endsToResize]) + IRanges::end(y)
        IRanges::start(tiles[endsToResize]) <- max(trimmedStarts, IRanges::start(y))
        IRanges::end(tiles[endsToResize]) <- IRanges::end(y)

        

        ## remove duplicate tiles if present
        tiles <- unique(tiles)
        return(tiles)
    }

    ## run lambda function
    tileList <- BiocParallel::bplapply(l$chromosomes, lambdaFun, sl = l)

    ## concatenate results into one list
    tiles <- do.call("c", tileList)
    GenomeInfoDb::seqlengths(tiles) <- GenomeInfoDb::seqlengths(l$chromosomes)
    GenomeInfoDb::seqlevels(tiles, force = TRUE) <- GenomeInfoDb::seqlevelsInUse(tiles)

    ## resize start and end tiles till correct tile size is reached
    startsToResize <- which(width(tiles) < l$tileSize & start(tiles) == 1)
    tiles[startsToResize] <- resize(tiles[startsToResize], width = l$tileSize)
    endsToResize <- which(width(tiles) < l$tileSize)
    tiles[endsToResize] <- resize(tiles[endsToResize], width = l$tileSize, fix = "end")
    
    ## add 'id' column, check element and put settings in metadata
    S4Vectors::mcols(tiles)$id <- 1:length(tiles)
    l$numTiles <- length(tiles)

    futile.logger::flog.trace(paste0("GenoGAMDataSet: Data split into ", l$numTiles, " tiles"))

    l$check <- TRUE
    S4Vectors::metadata(tiles) <- l
    return(tiles)
}


#' Convert the config columns to the right type.
#'
#' @param config A data.frame with pre-specified columns.
#' @param directory The directory of the files
#' @return The same data.frame with the columns of the right type.
#' @noRd
.normalizeConfig <- function(config, directory) {
    
    if(class(config) == "character") {
        config <- fread(config, header = TRUE, data.table = FALSE)
    }

    config$ID <- as.factor(config$ID)
    config$file <- file.path(directory, as.character(config$file))
    config$paired <- as.logical(config$paired)
   
    return(config)
}

## #' Construct GenomicTiles from a config file or a config data.frame
## #'
## #' See GenomicTiles in GenomicTiles-class.R for description.
## #' @return A GenomicTiles object.
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## .GenoGAMDataSetFromConfig <- function(config, chunkSize, overhangSize,
##                                     design, directory, settings, hdf5, ...) {

##     ## initialize some variables
##     args <- list()
        
##     ## normalize config object
##     config <- .normalizeConfig(config, directory)

##     center <- slot(settings, "center")
##     if(!is.null(center)) {
##         settings <- slot(settings, "processFunction") <- .processCountChunks
##     }

##     ## get chromosomeLengths
##     header <- Rsamtools::scanBamHeader(config$file[1])
##     chroms <- header[[1]]$targets
##     chromosomeLengths <- GenomicRanges::GRanges(names(chroms), IRanges(start = rep(1, length(chroms)),
##                                                         end = chroms))
##     GenomeInfoDb::seqlengths(chromosomeLengths) <- chroms
##     chromosomeList <- getDefaults(settings, "chromosomeList")
##     if(!is.null(chromosomeList)) {
##         GenomeInfoDb::seqlevels(chromosomeLengths, force = TRUE) <- chromosomeList
##         chroms <- chroms[names(chroms) %in% chromosomeList]
##     }

##     ## read in data
##     futile.logger::flog.info("Reading in data.")
##     rawData <- lapply(1:nrow(config), function(ii) {
##         if(!is.null(center)) {
##             args <- list(paired = config$paired[ii],
##                          center = center)
##         }
##         unlist(do.call(.readData,
##                       c(list(path = config$file[ii],
##                              processFUN =
##                              getDefaults(settings, "processFunction"),
##                              chromosomeList = chromosomeList,
##                              params = getDefaults(settings, "bamParams"),
##                              asMates = config$paired[ii]),
##                         args, list(...))), use.names = FALSE)
##     })

##     names(rawData) <- config$ID
##     assays <- DataFrame(rawData)

##     ## generate rowRanges
##     bamParamsWhich <- Rsamtools::bamWhich(getDefaults(settings, "bamParams"))
##     if(length(bamParamsWhich) != 0) {
##         gp <- GenomicRanges::GPos(GenomicRanges::GRanges(bamParamsWhich))
##         GenomeInfoDb::seqlengths(gp) <- chroms[GenomeInfoDb::seqlevels(GRanges(bamParamsWhich))]
##     }
##     else {
##         gp <- GenomicRanges::GPos(chromosomeLengths)
##         GenomeInfoDb::seqlengths(gp) <- chroms
##     }
    
##     gt <- GenomicTiles(assays, rowRanges = gp, chunkSize = chunkSize,
##                        overhangSize = overhangSize)
##     metadata(slot(gt, "index"))$check <- TRUE

##     ## make colData
##     colData <- DataFrame(config)[,-c(1:3), drop = FALSE]
##     rownames(colData) <- config$ID
##     colData(gt) <- colData
##     ## initiate size factors
##     sf <- rep(0, ncol(gt))
##     names(sf) <- colnames(gt)

##     settings <- setDefaults(settings, chromosomeList = GenomeInfoDb::seqlevels(tileSettings(gt)$chromosomes))

##     gtd <- new("GenoGAMDataSet", gt, settings = settings,
##                design = design, sizeFactors = sf)

##     ## check if everything was set fine
##     correct <- checkSettings(gtd)
##     if(!correct) break
    
##     futile.logger::flog.info("DONE")
##     return(gtd)
## }

## #' Make an example /code{GenoGAMDataSet}
## #'
## #' @return A /code{GenoGAMDataSet} object
## #' @examples
## #' test <- makeTestGenoGAMDataSet()
## #' @export

## Check functions
## ===============

#' Check a specified setting of GenoGAMDataSet
#'
#' @noRd
.checkSettings <- function(object, params = c("chunkSize", "tileSize",
                                       "equality", "tileRanges",
                                       "chromosomes", "numTiles")) {
    param <- match.arg(params)
    switch(param,
           chunkSize = .checkChunkSize(object),
           tileSize = .checkTileSize(object),
           equality = .checkEqualityOfTiles(object),
           chromosomes = .checkChromosomes(object),
           numTiles = .checkNumberOfTiles(object),
           tileRanges = .checkTileRanges(object))
}

## check the chunk size and return a logical value
.checkChunkSize <- function(object) {
    widths <- IRanges::width(getIndex(object))
    diffs <- (widths - 2*getOverhangSize(object)) - getChunkSize(object)
    res <- all.equal(diffs, rep(0, length(diffs)), tolerance = 0)
    return(res)
}

## check the tile size and return a logical value
.checkTileSize <- function(object) {
    widths <- IRanges::width(getIndex(object))
    diffs <- widths - getTileSize(object)
    res <- all.equal(diffs, rep(0, length(diffs)), tolerance = 0)
    return(res)
}

## check the overhang size and return a logical value
.checkEqualityOfTiles <- function(object) {
    widths <- width(getIndex(object))
    res <- all.equal(min(widths), max(widths))
    return(res)
}

## check the chromosome list and return a logical value
.checkChromosomes <- function(object) {
    objChroms <- GenomeInfoDb::seqlengths(object)
    indexChroms <- GenomeInfoDb::seqlengths(getIndex(object))
    validChroms <- GenomeInfoDb::seqlengths(getChromosomes(object))
    
    objChroms <- objChroms[order(names(objChroms))]
    indexChroms <- indexChroms[order(names(indexChroms))]
    validChroms <- validChroms[order(names(validChroms))]
    
    res1 <- all.equal(indexChroms, validChroms, objChroms)
    res2 <- all.equal(names(indexChroms), names(validChroms), names(objChroms))
    if(is(res1, "character") | is(res2, "character")) {
        ans <- paste("Chromosome Lengths:", res1, "\n", "Chromosome Names:", res2, "\n")
        return(ans)
    }
    return(res1 & res2)
}

## check the number of tiles and return a logical value
.checkNumberOfTiles <- function(object) {
    tiles <- length(getIndex(object))
    validTiles <- getTileNumber(object)
    res <- all.equal(tiles, validTiles)
    return(res)
}

## check ranges
.checkTileRanges <- function(object) {
    tileRanges <- tileSettings(object)$chromosomes
    dataRanges <- dataRange(object)
    res <- all.equal(tileRanges, dataRanges)
    return(res)
}

##' Function to check the GenoGAMDataSet object
##'
##' @noRd
.checkGenoGAMDataSet <- function(object) {
    futile.logger::flog.trace("Check if tile settings match the data.")

    params = c("chunkSize", "tileSize", "tileRanges",
               "equality", "chromosomes", "numTiles")

    settings <- tileSettings(object)
    
    if(is.null(settings$check)) {
        futile.logger::flog.warn("Checks dismissed due to empty object or forgotten setting")
        return(FALSE)
    }
    if(!settings$check) {
        futile.logger::flog.warn("Settings checking deactivated. Modeling on these tiles might yield wrong results.")
        return(FALSE)
    }
    res <- sapply(params, .checkSettings, object = object)
    if(is(res, "character")) {
        errorIndx <- which(res != "TRUE")
        futile.logger::flog.error("Checks failed. Following settings display errors:\n")
        print(res[errorIndx])
        return(FALSE)
    }
    
    futile.logger::flog.trace("All checks passed.")
    return(TRUE)
}

#' @noRd
setGeneric("checkObject", function(object) standardGeneric("checkObject")) 


#' Check data compliance with tile settings
#'
#' Check if the indices were build correctly, according to the
#' specified parameters
#'
#' @rdname checkObject
#' @param object A /code{GenomicTiles} object.
#' @return A logical value
#' @examples 
#' gt <- makeTestGenomicTiles()
#' checkSettings(gt)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
setMethod("checkObject", "GenoGAMDataSet", function(object) {
    .checkGenoGAMDataSet(object)
})


## Accessors
## =========

##' @export
setGeneric("getIndex", function(object) standardGeneric("getIndex"))


##' @describeIn GenoGAMDataSet An accessor to the index slot
setMethod("getIndex", signature(object = "GenoGAMDataSet"), function(object) {
    return(slot(object, "index"))
})


##' @export
setGeneric("tileSettings", function(object) standardGeneric("tileSettings"))

##' @describeIn GenoGAMDataSet The accessor to the list of settings,
##' that were used to generate the tiles.
setMethod("tileSettings", "GenoGAMDataSet", function(object) {
    metadata(getIndex(object))
})

##' @export
setGeneric("dataRange", function(object) standardGeneric("dataRange")) 

##' @describeIn GenoGAMDataSet The actual underlying GRanges showing
##' the range of the data.
setMethod("dataRange", "GenoGAMDataSet", function(object) {
    rowRanges(object)@pos_runs
})

##' @export 
setGeneric("getChromosomes", function(object) standardGeneric("getChromosomes")) 

##' @describeIn GenoGAMDataSet A GRanges object representing the chromosomes
##' or chromosome regions on which the model will be computed
setMethod("getChromosomes", "GenoGAMDataSet", function(object) {
    tileSettings(object)$chromosomes
})

##' @export 
setGeneric("getTileSize", function(object) standardGeneric("getTileSize")) 

##' @describeIn GenoGAMDataSet The size of the tiles
setMethod("getTileSize", "GenoGAMDataSet", function(object) {
    tileSettings(object)$tileSize
})

##' @export 
setGeneric("getChunkSize", function(object) standardGeneric("getChunkSize")) 

##' @describeIn GenoGAMDataSet The size of the chunks
setMethod("getChunkSize", "GenoGAMDataSet", function(object) {
    tileSettings(object)$chunkSize
})

##' @export 
setGeneric("getOverhangSize", function(object) standardGeneric("getOverhangSize")) 

##' @describeIn GenoGAMDataSet The size of the overhang (on one side)
setMethod("getOverhangSize", "GenoGAMDataSet", function(object) {
    tileSettings(object)$overhangSize
})

##' @export 
setGeneric("getTileNumber", function(object) standardGeneric("getTileNumber")) 

##' @describeIn GenoGAMDataSet The total number of tiles
setMethod("getTileNumber", "GenoGAMDataSet", function(object) {
    tileSettings(object)$numTiles
})

##' @describeIn GenoGAMDataSet Access to the design slot.
setMethod("design", "GenoGAMDataSet", function(object) {
    slot(object, "design")
})

##' @describeIn GenoGAMDataSet Replace method of the design slot.
setReplaceMethod("design", "GenoGAMDataSet", function(object, value) {
    slot(object, "design") <- value
    return(object)
})

##' @describeIn GenoGAMDataSet Access to the sizeFactors slot
setMethod("sizeFactors", "GenoGAMDataSet", function(object) {
    slot(object, "sizeFactors")
})

##' @describeIn GenoGAMDataSet Replace method of the sizeFactors slot
setReplaceMethod("sizeFactors", "GenoGAMDataSet", function(object, value) {
    slot(object, "sizeFactors") <- value
    return(object)
})












## makeTestGenoGAMDataSet <- function() {
##     gp <- GenomicRanges::GPos(GenomicRanges::GRanges(c("chrI", "chrII"), IRanges(c(1,1), c(50,50))))
##     df <- DataFrame(a = Rle(1:100), b = Rle(101:200))
##     se <- SummarizedExperiment::SummarizedExperiment(list(df), rowRanges = gp)
##     ggd <- GenoGAMDataSet(se, chunkSize = 15, overhangSize = 3,
##                           design = ~s(x))
## }

## ## Coercion
## ## ========

## # converts the GenomicTiles object to a DataFrame
## .GenoGAMDataSetToDataFrame <- function(from) {
##     df <- as.data.frame(rowRanges(from))
##     numDataFrames <- length(assays(from))
##     gtDim <- dim(from)
    
##     res <- unlist(assays(from))   
##     gtdf <- DataFrame(cbind(df, res))

##     metadata(gtdf) <- list(sizeFactors = sizeFactors(from))
        
##     return(gtdf)
## }

## #' GenoGAMDataSet to DataFrame
## #' 
## #' @name GenoGAMDataSetToDataFrame
## #' @family GenoGAMDataSet
## setAs(from = "GenoGAMDataSet", to = "DataFrame", def = function(from) {
##     .GenoGAMDataSetToDataFrame(from)
## })

## setMethod("as.data.frame", "GenoGAMDataSet", function(x) {
##     ## Bug in as(x, "data.frame"), replace with
##     as.data.frame(as(x, "DataFrame"))
## })

## ## Subsetting
## ## ==========
## #' Subset method for \code{GenoGAMDataSet}
## #'
## #' Subsetting the \code{GenoGAMDataSet} by a logical statement
## #'
## #' @param x A \code{GenoGAMDataSet} object.
## #' @param ... Further arguments. Mostly a logical statement.
## #' Note that the columnnames for chromosomes and positions
## #' are: seqnames and pos.
## #' @return A subsetted \code{GenomicTiles} object.
## #' @examples
## #' ggd <- makeTestGenoGAMDataSet()
## #' res <- subset(ggd, seqnames == "chrI" & pos <= 50)
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## setMethod("subset", "GenoGAMDataSet", function(x, ...) {
##     settings <- getSettings(x)
##     design <- design(x)
##     sf <- sizeFactors(x)
##     gt <- .subsetGenomicTiles(x, ...)
##     settings <- setDefaults(settings, chromosomeList = GenomeInfoDb::seqlevels(gt))

##     gtd <- new("GenoGAMDataSet", gt, settings = settings,
##                design = design, sizeFactors = sf)
##     return(gtd)
## })

## #' Subset by overlaps method for \code{GenoGAMDataSet}
## #'
## #' Subsetting the \code{GenoGAMDataSet} by a \code{GRanges} object
## #'
## #' @param query A \code{GenoGAMDataSet} object.
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
## #' @return A subsetted \code{GenoGAMDataSet} object.
## #' @examples
## #' ggd <- makeTestGenoGAMDataSet()
## #' gr <- GRanges("chrI", IRanges(1,50))
## #' res <- subsetByOverlaps(ggd, gr)
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## setMethod("subsetByOverlaps", c("GenoGAMDataSet", "GRanges"),
##           function(query, subject, maxgap=0L, minoverlap=1L,
##                       type=c("any", "start", "end", "within", "equal"),...) {
##               settings <- getSettings(query)
##               design <- design(query)
##               sf <- sizeFactors(query)
##               subgt <- .subsetByOverlaps(query, subject, maxgap=maxgap,
##                                          minoverlap=minoverlap,
##                                          type=type,...)
##               settings <- setDefaults(settings, chromosomeList = GenomeInfoDb::seqlevels(subgt))
              
##               gtd <- new("GenoGAMDataSet", subgt, settings = settings,
##                          design = design, sizeFactors = sf)
##               return(gtd)
##           })

## #' Subsetting by GRanges
## #'
## #' Providing subsetting by GRanges through the single-bracket operator
## #'
## #' @param x A \code{GenoGAMDataSet} object
## #' @param i A \code{GRanges} object
## #' @return A subsetted \code{GenoGAMDataSet} object
## #' @rdname GenoGAMDataSet-brackets
## setMethod("[", c("GenoGAMDataSet", "GRanges"), function(x, i) {
##     settings <- getSettings(x)
##     design <- design(x)
##     sf <- sizeFactors(x)
##     subgt <- .exactSubsetByOverlaps(x, i)
##     settings <- setDefaults(settings, chromosomeList = GenomeInfoDb::seqlevels(subgt))
              
##     gtd <- new("GenoGAMDataSet", subgt, settings = settings,
##                design = design, sizeFactors = sf)
##     return(gtd)
## })
    

## ## Silent methods
## ## ==============
## setGeneric("getSettings", function(object) standardGeneric("getSettings"))

## setMethod("getSettings", "GenoGAMDataSet", function(object) {
##     slot(object, "settings")
## })
