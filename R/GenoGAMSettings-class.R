## ===============
## GenoGAMSettings
## ===============

setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("logicalOrNULL", c("logical", "NULL"))
setClassUnion("functionOrNULL", c("function", "NULL"))

#' GenoGAMSettings
#'
#' This class is designed to store global settings for the computation of the
#' GenoGAM package
#'
#' @param ... Any parameters corresponding to the slots and their possible
#' values.
#' 
#' @slot center A logical or NULL value to specify if the raw data should be
#' centered, i.e. only the midpoint of the fragment will be used to
#' represent its coverage. See details.
#' @slot chromosomeList A character vector of chromosomes to be used.
#' NULL for all chromosomes.
#' @slot bamParams An object of class ScanBamParam.
#' See ?Rsamtools::ScanBamParam for possible settings. Usually used to set
#' specific ranges, to read in.
#' @slot processFunction A custom function on how to process raw data. Not
#' used if center is TRUE/FALSE. This is not intended for the user, but if
#' needed anyway, see details.
#' @slot optimMethod The optiomisation method to be used in cross validation.
#' See ?optim for a complete list.
#' @slot optimControl Settings for the optim function. See ?optim for a 
#' complete list.
#' @details Center can have three values: TRUE, FALSE, NULL. TRUE will
#' trigger the center function, FALSE will trigger the use of the entire
#' fragment. NULL should be used in case a custom process function is used.
#' In case a custom function is used, it has to satisfy the following:
#' It has to handle a GAlignments object as input and output a 
#' GRanges object of regions, e.g. fragments. This regions are in turn used
#' to compute the coverage via the IRanges::coverage function. Note, that there
#' is a difference between the GAlignments object in the single and paired end
#' case.
#' @name GenoGAMSettings-class
#' @rdname GenoGAMSettings-class
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @exportClass GenoGAMSettings
setClass("GenoGAMSettings",
         slots = list(center = "logicalOrNULL", chromosomeList = "characterOrNULL",
             bamParams = "ScanBamParam", processFunction = "functionOrNULL", 
             optimMethod = "character", optimControl = "list"),
         prototype = list(center = TRUE, chromosomeList = NULL,
             bamParams = Rsamtools::ScanBamParam(what = c("pos", "qwidth")),
             processFunction = NULL, optimMethod = "Nelder-Mead",
             optimControl = list(maxit = 50, fnscale = -1, trace = 10)))

## Validity
## ========

#' Validating the correct type

.validateBAMParamsType <- function(object) {
    if(class(object@bamParams) != "ScanBamParam") {
        return("'bamParams' must be a ScanBamParam object")
    }
    NULL
}

.validateOptimMethodType <- function(object) {
    if(class(object@optimMethod) != "character") {
        return("'optimMethod' must be a character object")
    }
    NULL
}

.validateOptimControlType <- function(object) {
    if(class(object@optimControl) != "list") {
        return("'optimControl' must be a list object")
    }
    NULL
}

## general validate function
.validateGenoGAMSettings <- function(object) {
    c(.validateBAMParamsType(object), .validateOptimMethodType(object),
      .validateOptimControlType(object))
}

setValidity2("GenoGAMSettings", .validateGenoGAMSettings)

## Constructor
## ==========

#' The Constructor function for GenoGAMSettings is merely a wrapper for
#' new("GenoGAMSettings", ...)
#' 
#' @name GenoGAMSettings
#' @rdname GenoGAMSettings-class
#' @export
GenoGAMSettings <- function(...) {
    return(new("GenoGAMSettings", ...)) 
}

## Cosmetics
## =========

## show method for GenoGAMSettings
.showGenoGAMSettings <- function(ggs) {
    if(!is.null(ggs@center)) {
        processFun <- FALSE
    }
    else {
        processFun <- TRUE
    }
    
    cat("-------------------- Read-in parameters -----------------\n")
    cat("Center:", ggs@center, "\n")
    cat("Chromosomes:", paste(ggs@chromosomeList, collapse = ", "), "\n")
    cat("Custom process function:", ifelse(processFun, "enabled,", "disabled"), 
        "\n")
    cat("\n")
    cat("-------------------- BAM parameters ---------------------\n")
    show(ggs@bamParams)
    cat("\n")
    cat("-------------------- Parallel backend -------------------\n")
    show(BiocParallel::registered()[[1]])
    cat("\n")
    cat("-------------------- Cross Validation parameters --------\n")
    cat("Optimization method:", ggs@optimMethod, "\n")
    cat("Optimization control:\n")
    for(ii in 1:length(ggs@optimControl)) {
        cat(paste0("  ", names(ggs@optimControl)[ii], ": ",
                   ggs@optimControl[[ii]], "\n"))
    }
}

setMethod("show", "GenoGAMSettings", function(object) {
    .showGenoGAMSettings(object)
})
