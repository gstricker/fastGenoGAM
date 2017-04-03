#############################
## GenoGAM class 
##################

#' @include GenoGAMSettings-class.R
NULL

#' GenoGAM class
#'
#' This is the class the holds the complete model as well all hyperparameters
#' and settings that were used to fit it. It extends the RangedSummarizedExperiment
#' class by adding a couple of more slots to hold hyperparameters and settings.
#' The 'assays' slot holds the basepair fit and standard deviation. Additionally
#' all knot positions and beta coefficients will be stored in the 'smooths' slot
#' in order to be able to make use of the piecewise function that produces the
#' fit. For information on the slots inherited from SummarizedExperiment
#' check the respective class.
#' 
#' @slot family The name of the distribution family used
#' @slot design The formula of the model
#' @slot sizeFactors The offset used in the model. 
#' @slot factorialDesign The factorial design used. The same as colData in the
#' GenoGAMDataSet
#' @slot params All hyperparameters used to fit the data. The parameters
#' estimated by cross validation can also be found here. But the parameters
#' used in cross validation are in the settings slot.
#' @slot settings A GenoGAMSettings object representing the global
#' settings that were used to compute the model.
#' @slot smooths A data.table of knots and coefficients of the model
#' @name GenoGAM-class
#' @rdname GenoGAM-class
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @exportClass GenoGAM
setClass("GenoGAM",
         contains = "RangedSummarizedExperiment",
         slots = list(family = "character",
                      design = "formula",
                      sizeFactors = "numeric",
                      factorialDesign = "DataFrame",
                      params = "list",
                      settings = "GenoGAMSettings",
                      smooths = "data.table"),
         prototype = prototype(family = "nb",
                               design = ~ s(x),
                               sizeFactors = numeric(),
                               factorialDesign = DataFrame(),
                               params = list(),
                               settings = GenoGAMSettings(),
                               smooths = data.table::data.table()))

## Validity
## ========

.validateFamilyType <- function(object) {
    if(class(slot(object, "family")) != "character") {
        return("'family' must be a character object")
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

.validateFactorialDesign <- function(object) {
    if(class(slot(object, "factorialDesign")) != "DataFrame") {
        return("'factorialDesign' must be a DataFrame class")
    }
}

.validateParamsType <- function(object) {
    if(class(slot(object, "params")) != "list") {
        return("'params' must be a list object")
    }
    NULL
}

.validateSettingsType <- function(object) {
    if(class(slot(object, "settings")) != "GenoGAMSettings") {
        return("'settings' must be a GenoGAMSettings object")
    }
    NULL
}

.validateSmoothsType <- function(object) {
    if(class(slot(object, "smooths"))[1] != "data.table") {
        return("'smooths' must be a data.table object")
    }
    NULL
}

## general validate function
.validateGenoGAM <- function(object) {
    c(.validateFamilyType(object),
      .validateDesignType(object),
      .validateSFType(object),
      .validateFactorialDesign(object),
      .validateParamsType(object),
      .validateSettingsType(object),
      .validateSmoothsType(object))
}

setValidity2("GenoGAM", .validateGenoGAM)


## Constructor
## ========

#' GenoGAM constructor
#'
#' The GenoGAM constructor, not designed to be actually used, by the user.
#' Rather to be a point of reference and documentation for slots and how
#' to access them.
#'
#' @aliases design sizeFactors getSettings getFamily colData getParams getSmooths
#' @param object,x For use of S4 methods. The GenoGAM object.
#' @param ... Slots of the GenoGAM class. See the slot description.
#' @return An object of the type GenoGAM.
#' @name GenoGAM
#' @rdname GenoGAM-class
#' @export
GenoGAM <- function(...) {
    return(new("GenoGAM", ...))
}


## Accessors
## =========

##' @describeIn GenoGAM An accessor to the design slot
setMethod("design", "GenoGAM", function(object) {
    slot(object, "design")
})

##' @describeIn GenoGAM An accessor to the sizeFactors slot
setMethod("sizeFactors", "GenoGAM", function(object) {
    slot(object, "sizeFactors")
})

##' @export
setGeneric("getSettings", function(object) standardGeneric("getSettings"))

##' @describeIn GenoGAM An accessor to the settings slot
setMethod("getSettings", "GenoGAM", function(object) {
    slot(object, "settings")
})

##' @export
setGeneric("getFamily", function(object) standardGeneric("getFamily"))

##' @describeIn GenoGAM An accessor to the family slot
setMethod("getFamily", "GenoGAM", function(object) {
    slot(object, "family")
})

##' @describeIn GenoGAM An accessor to the factorialDesign slot.
##' It overwrites the inherited colData function, which points to the colData slot.
setMethod("colData", "GenoGAM", function(x) {
    slot(x, "factorialDesign")
})

##' @export
setGeneric("getParams", function(object) standardGeneric("getParams"))

##' @describeIn GenoGAM An accessor to the params slot
setMethod("getParams", "GenoGAM", function(object) {
    slot(object, "params")
})

##' @export
setGeneric("getSmooths", function(object) standardGeneric("getSmooths"))

##' @describeIn GenoGAM An accessor to the smooths slot
setMethod("getSmooths", "GenoGAM", function(object) {
    slot(object, "smooths")
})


## ## Cosmetics
## ## ==========

## #' The actual show function
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## .showGenoGAM <- function(gg) {
##     fitparams <- slot(gg, "fitparams")
##     cvparams <- slot(gg, "cvparams")
##     settings <- metadata(getIndex(gg))

##     cat("Family: ")
##     show(slot(gg, "family"))
##     cat("Formula:\n")
##     design(gg)
##     cat("\n")
    
##     cat("Experiment Design:\n")
##     colData(gg)
##     cat("\n")
##     cat("Size factors: \n")
##     sizeFactors(gg)
##     cat("\n")
##     cat("Global Parameters:\n")
##     cat("  Lambda: ", fitparams["lambda"], "\n")
##     cat("  Theta: ", fitparams["theta"], "\n")
##     cat("  Coefficient of Variation: ", fitparams["CoV"], "\n")
##     cat("  H: ", fitparams["H"], "\n")
##     cat("  Effective degrees of freedom: ", fitparams["df"], "\n")
##     cat("  Approximate windown size: ", settings$tileSize/fitparams["df"], "\n")

##     cat("\n")
##     cat("Cross Validation:")
##     if(!is.na(cvparams["cv"]) & cvparams["cv"] > 0) {
##         cat(" Performed\n")
##     }
##     else {
##         cat(" Not performed\n")
##     }
##     cat("  K-folds:", cvparams["kfolds"], "\n")
##     cat("  Number of tiles:", cvparams["ncv"], "\n")
##     cat("  Interval size:", cvparams["size"], "\n")

##     cat("\n")
##     cat("Tile settings:\n")
##     cat("  chunk size:", settings$chunkSize, "\n")
##     cat("  tile size:", settings$tileSize, "\n")
##     cat("  overhang size:", settings$overhangSize, "\n")
##     cat("  number of tiles:", settings$numTiles, "\n")
## }

## setMethod("show", "GenoGAM", function(object) {
##     .showGenoGAM(object)
## })


## setMethod("summary", "GenoGAM", function(object) {
##     .showGenoGAM(object)
##     cat("Coefficients: \n")
##     coef(object)
## })

## #' Make an example /code{GenoGAM}
## #'
## #' @return A /code{GenoGAM} object
## #' @examples
## #' test <- makeTestGenoGAM()
## #' @export
## makeTestGenoGAM <- function() {
##     gr <- GenomicRanges::GRanges("chrI", IRanges::IRanges(1, 100))
##     gp <- GenomicRanges::GPos(gr)
##     fits <- DataFrame(runif(100), runif(100), runif(100), runif(100))
##     names(fits) <- c("s(x)", "s(x):type", "se.s(x)", "se.s(x):type")
##     se <- SummarizedExperiment::SummarizedExperiment(rowRanges = gp, 
##                                                      assays = list(fits))
##     gt <- GenomicTiles(se, chunkSize = 10) 
##     fitparams <- c(lambda = 10, theta = 2, CoV = 0.2, H = 0.001, df = 3.6)
##     cvparams <- c(cv = FALSE, kfolds = 10, ncv = 20, size = 20)
##     gg <- new("GenoGAM", gt)
##     slot(gg, "fitparams") <- fitparams
##     slot(gg, "cvparams") <- cvparams
##     return(gg)
## }












## ## Cosmetics
## ## ==========

## #' The actual show function
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## .showGenoGAM <- function(gg) {
##     fitparams <- slot(gg, "fitparams")
##     cvparams <- slot(gg, "cvparams")
##     settings <- slot(gg, "tileSettings")
    
##     show(slot(gg, "family"))
##     cat("Formula:\n")
##     show(slot(gg, "design"))
##     cat("\n")
    
##     cat("Experiment Design:\n")
##     show(slot(gg, "experimentDesign"))
##     cat("\n")
##     cat("Global Estimates:\n")
##     cat("  Lambda:", fitparams["lambda"], "\n")
##     cat("  Theta:", fitparams["theta"], "\n")
##     cat("  Coefficient of Variation:", fitparams["CoV"], "\n")

##     cat("\n")
##     cat("Cross Validation:")
##     if(!is.na(cvparams["cv"]) & cvparams["cv"] > 0) cat(" Performed\n")
##     else cat(" Not performed\n")
##     cat("  K-folds:", cvparams["kfolds"], "\n")
##     cat("  Number of tiles:", cvparams["ncv"], "\n")
##     cat("  Interval size:", cvparams["size"], "\n")

##     cat("\n")
##     cat("Tile settings:\n")
##     cat("  chunk size:", settings$chunkSize, "\n")
##     cat("  tile size:", settings$tileSize, "\n")
##     cat("  overhang size:", settings$overhangSize, "\n")
##     cat("  number of tiles:", settings$numTiles, "\n")
## }

## setMethod("show", "GenoGAM", function(object) {
##     .showGenoGAM(object)
## })

## setMethod("summary", "GenoGAM", function(object) {
##     fits <- data.table::data.table(cbind(seqnames(slot(object, "positions")),
##                   GenomicRanges::pos(slot(object, "positions")),
##                   slot(object, "fits")))
##     data.table::setnames(fits, names(fits)[1:2], c("seqnames", "x"))
##     .showGenoGAM(object)
##     cat("\n")
##     cat("Fitted values:\n")
##     cat("\n")
##     show(fits)
## })

## #' Make an example /code{GenoGAM}
## #'
## #' @return A /code{GenoGAM} object
## #' @examples
## #' test <- makeTestGenoGAM()
## #' @export
## makeTestGenoGAM <- function() {
##     gg <- .GenoGAM()
##     gp <- GenomicRanges::GPos(GenomicRanges::GRanges("chrI", IRanges::IRanges(1, 100)))
##     fits <- data.frame(runif(100), runif(100), runif(100), runif(100))
##     names(fits) <- c("s(x)", "s(x):type", "se.s(x)", "se.s(x):type")
##     slot(gg, "positions") <- gp
##     slot(gg, "fits") <- fits
##     return(gg)
## }

## .subsetByRanges <- function(gg, ranges) {
##     pos <- rowRanges(gg)
##     ov <- findOverlaps(pos, ranges)
##     indx <- queryHits(ov)
##     slot(gg, "positions") <- pos[indx,]
##     slot(gg, "fits") <- slot(gg, "fits")[indx,]
##     return(gg)
## }

## .subsetByPosition <- function(gg, ...) {
##     positions <- slot(gg, "positions")
##     ov <- findOverlaps(positions, subset(positions, ...))
##     indx <- queryHits(ov)
##     slot(gg, "positions") <- positions[indx,]
##     slot(gg, "fits") <- slot(gg, "fits")[indx,]
##     return(gg)
## }
## #' Subset method for \code{GenoGAM}
## #'
## #' Subsetting the \code{GenoGAM} by a logical statement
## #'
## #' @param x A \code{GenoGAM} object.
## #' @param ... Further arguments. Mostly a logical statement.
## #' Note that the columnnames for chromosomes and positions
## #' are: seqnames and pos.
## #' @return A subsetted \code{GenoGAM} object.
## #' @examples
## #' gg <- makeTestGenoGAM()
## #' subset(gg, pos <= 40)
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## setMethod("subset", "GenoGAM", function(x, ...) {
##     .subsetByPosition(x, ...)
## })

## #' Subset by overlaps method for \code{GenoGAM}
## #'
## #' Subsetting the \code{GenoGAM} by a \code{GRanges} object
## #'
## #' @param query A \code{GenoGAM} object.
## #' @param subject A \code{GRanges} object
## #' @param ... Additional parameters
## #' @return A subsetted \code{GenoGAM} object.
## #' @examples
## #' gg <- makeTestGenoGAM()
## #' gr <- GRanges("chrI", IRanges(1,40))
## #' subsetByOverlaps(gg, gr)
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## setMethod("subsetByOverlaps", "GenoGAM", function(query, subject) {
##     .subsetByRanges(query, subject)
## })

## #' View the dataset
## #'
## #' Cbinding the columns all together and coercing to data.frame
## #'
## #' @param object A \code{GenoGAM} object
## #' @param ranges A \code{GRanges} object. Makes it possible to
## #' select regions by \code{GRanges}. Either ranges or seqnames, start and
## #' end must be supplied
## #' @param seqnames A chromosomes name. Either ranges or seqnames, start and
## #' end must be supplied
## #' @param start A start site. Either ranges or seqnames, start and
## #' end must be supplied
## #' @param end An end site. Either ranges or seqnames, start and
## #' end must be supplied
## #' @return A data.frame of the selected data.
## #' @examples
## #' gg <- makeTestGenoGAM()
## #' gr <- GRanges("chrI", IRanges(1,40))
## #' head(view(gg, gr))
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## #' @rdname GenoGAM-view
## #' @export
## setMethod("view", "GenoGAM", function(object, ranges = NULL, seqnames = NULL,
##                                       start = NULL, end = NULL) {
##     ## keep "seqnames" for consistency with Bioc, but rename variable as subset in
##     ## .subsetByPosition does not work properly otherwise
##     chromosome <- seqnames 
##     if(is.null(seqnames) & is.null(start) & is.null(end) & is.null(ranges)) {
##         temp <- object
##     }
##     else {
##         if(!is.null(ranges)) {
##             temp <- .subsetByRanges(object, ranges)
##         }
##         else {
##             if(is.null(start)) {
##                 start <- 1
##             }
##             if(is.null(end)) {
##                 end <- Inf
##             }
##             if(is.null(seqnames)){
##                 temp <- .subsetByPosition(object, pos >= start & pos <= end)
##             }
##             else {
##                 temp <- .subsetByPosition(object, seqnames == chromosome & pos >= start & pos <= end)
##             }
##         }
##     }
##     res <- cbind(slot(temp, "positions"), slot(temp, "fits"))
##     return(res)
## })

## .pvals <- function(gg, log.p = FALSE) {
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

## .result <- function(gg, log.p = FALSE) {
##     if(nrow(getFits(gg)) > 0) {
##         pvals <- .pvals(gg, log.p)
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

## plot.GenoGAM <- function(object, ranges = NULL, seqnames = NULL,
##                          start = NULL, end = NULL, scale = TRUE,
##                          pages = 1, select = NULL) {
##     base <- TRUE
##     if(require(ggplot2)) {
##         base <- FALSE
##     }
##     sub <- view(ranges = ranges, seqnames = seqnames,
##                 start = start, end = end)
##     if(is.null(select)) {
##         select <- c("s(x)", paste("s(x)", colnames(design(object)), sep = ":"))
##     }
##     else {
##         select <- paste("s(x)", select, sep = ":")
##     }

##     if(base) {
##         plot_base(sub, scale = scale, pages = pages,
##                   select = select)
##     }
## }

## ## #' Prediction from fitted GenoGAM
## ## #'
## ## #' Takes a fitted ‘GenoGAM’ object produced by ‘genogam()’ and returns
## ## #' predictions given a new set of values for the model covariates or
## ## #' the original values used for the model fit. Predictions can be
## ## #' accompanied by standard errors, based on the posterior distribution
## ## #' of the model coefficients.
## ## #'
## ## #' @param object A GenoGAM object.
## ## #' @param newdata A GenoGAMDataSet object parallel to the original data set.
## ## #' @param type Either 'terms', 'link' or 'response'. The first returns the
## ## #' fits for each spline function. The last two return the fit for the response
## ## #' variable transformed by the link function ('link') or not ('response').
## ## #' @param se.fit Should the standard errors be returned?
## ## #' @param terms Which terms should be returned?
## ## #' @param exclude Which terms should be excluded?
## ## #' @return A list of fits and optinally of standard errors.
## ## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## ## #' @export
## ## predict.GenoGAM <- function(object, newdata, type = c("terms", "response", "link"),
## ##                             se.fit = FALSE,
## ##                             terms = NULL, exclude = NULL) {
## ##     type <- match.arg(type)
## ##     if(missing(newdata)) x <- slot(object, "positions")
## ##     else x <- rowRanges(newdata)
## ##     rows <- match(pos(x), pos(slot(object, "positions")))

## ##     secols <- which(unlist(regexec("se.fit", names(object@fits))) >= 1)

## ##     if(is.null(terms)) {
## ##         cols <- which(unlist(regexec("se.fit", names(object@fits))) != 1)
## ##     }
## ##     else {
## ##         cols <- sapply(terms, function(y) {
## ##             which(unlist(regexec(y, names(object@fits))) >= 1)
## ##         })
## ##         if(se.fit) secols <- cols[cols %in% secols]
## ##         cols <- cols[!(cols %in% secols)]
## ##     }

## ##     if(!is.null(exclude)) {
## ##       excols <- sapply(exclude, function(y) {
## ##           which(unlist(regexec(y, names(object@fits))) >= 1)
## ##       })
## ##       cols <- cols[!(cols %in% excols)]
## ##       secols <- secols[!(secols %in% excols)]
## ##     }
    
## ##     if(type == "response" | type == "link") {
## ##         cols <- which(unlist(regexec("se.fit", names(object@fits))) != 1)
## ##         fits <- rowSums(slot(object, "fits")[rows, cols])
## ##         if(se.fit) sefits <- rowSums(slot(object, "fits")[rows, secols])
## ##         if(type == "response") {
## ##             fits <- slot(object, "family")$linkinv(fits)
## ##             if(se.fit) sefits <- slot(object, "family")$linkinv(sefits)
## ##         }
## ##     }
## ##     else {
## ##         fits <- slot(object, "fits")[rows, cols]
## ##         if(se.fit) sefits <- slot(object, "fits")[rows, secols]
## ##     }
## ##     if(se.fit) return(list(fit = fits, se.fit = sefits))
## ##     return(list(fit = fits))
## ## }
