#################
## Setup class ##
#################

#' GenoGAMSEtup class
#'
#' A class to embody the setup for a GenoGAM fit
#' 
#' @slot params A list of all hyper-parameters which are either estimated
#' or fixed.
#' At the moment the smoothing parameter lambda and a second regularization
#' matrix H are provided.
#' @slot knots A list of knot positions on each chromosome.
#' @slot designMatrix The design matrix.
#' @slot beta The vector of coefficients to be estimated. Initialized.
#' @slot vcov The covariance matrix
#' @slot penaltyMatrix The penalty matrix S with penalization order r. 
#' By default r = 2.
#' @slot formula The formula of the model. Usually the same as the design of
#' the GenoGAMDataSet
#' @slot offset An offset of the samples
#' @slot family The distribution to be used. At the moment only "nb"
#' (Negative Binomial) is available.
#' @slot response The response vector
#' @slot fits The vector of fits
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
setClass("GenoGAMSetup",
         slots = list(params = "list", knots = "list",
                      designMatrix = "dgCMatrix", beta = "matrix",
                      vcov = "dgCMatrix", penaltyMatrix = "dgCMatrix",
                      formula = "formula", offset = "numeric", 
                      family = "character", response = "numeric",
                      fits = "numeric"),
         prototype = list(params = list(lambda = 0, theta = 0, H = 0,
                                        order = 2, penorder = 2),
                          knots = list(), designMatrix = new("dgCMatrix"),
                          beta = matrix(), vcov = new("dgCMatrix"),
                          penaltyMatrix = new("dgCMatrix"), formula = ~1,
                          offset = numeric(), family = "nb", 
                          response = numeric(), fits = numeric()))

## Validity
## ========

#' Validating the correct type
.validateParamsType <- function(object) {
    if(class(slot(object, "params")) != "list") {
        return("'params' must be a list object")
    }
    NULL
}

.validateParamsElements <- function(object) {
    if(!all(names(slot(object, "params")) %in% c("lambda", "theta", "H", "order", "penorder"))) {
        return("'params' must contain the elements 'lambda', 'theta', 'H', 'order' and 'penorder'")
    }
    NULL
}

.validateKnotsType <- function(object) {
    if(class(slot(object, "knots")) != "list") {
        return("'knots' must be a list object")
    }
    NULL
}

.validateDesignMatrixType <- function(object) {
    if(class(slot(object, "designMatrix")) != "dgCMatrix") {
        return("'designMatrix' must be a dgCMatrix object")
    }
    NULL
}

.validateBetaType <- function(object) {
    if(class(slot(object, "beta")) != "matrix") {
        return("'beta' must be a matrix object")
    }
    NULL
}

.validateCovarianceType <- function(object) {
    if(class(slot(object, "vcov")) != "dgCMatrix") {
        return("'vcov' must be a dgCMatrix object")
    }
    NULL
}

.validatePenaltyMatrixType <- function(object) {
    if(class(slot(object, "penaltyMatrix")) != "dgCMatrix") {
        return("'penaltyMatrix' must be a dgCMatrix object")
    }
    NULL
}

.validateFormulaType <- function(object) {
    if(class(slot(object, "formula")) != "formula") {
        return("'formula' must be a formula object")
    }
    NULL
}

.validateOffsetType <- function(object) {
    if(mode(slot(object, "offset")) != "numeric") {
        return("'offset' must be a numeric object")
    }
    NULL
}

.validateFamilyType <- function(object) {
    if(class(slot(object, "family")) != "character") {
        return("'family' must be a character object")
    }
    NULL
}

.validateResponseType <- function(object) {
    if(mode(slot(object, "response")) != "numeric") {
        return("'response' must be a numeric object")
    }
    NULL
}

.validateFitsType <- function(object) {
    if(mode(slot(object, "fits")) != "numeric") {
        return("'fits' must be a numeric object")
    }
    NULL
}

## general validate function
.validateGenoGAMSetup <- function(object) {
    c(.validateParamsType(object), .validateKnotsType(object),
      .validateDesignMatrixType(object), .validateBetaType(object),
      .validateCovarianceType(object), .validatePenaltyMatrixType(object),
      .validateFormulaType(object), .validateOffsetType(object),
      .validateFamilyType(object), .validateResponseType(object),
      .validateFitsType(object), .validateParamsElements(object))
}

setValidity2("GenoGAMSetup", .validateGenoGAMSetup)

#' Constructor
#' @noRd
GenoGAMSetup <- function(...) {
    ggs <- new("GenoGAMSetup", ...)
    params <- slot(ggs, "params")
    coreParams <- c("lambda", "theta", "H", "order", "penorder")
    allin <- coreParams  %in% names(params)
    if(!all(allin)) {
        coreValues <- list(lambda = 0, theta = 0, H = 0, order = 2, penorder = 2)
        for(elem in coreParams) {
            if(is.null(params[[elem]])) {
                params[[elem]] <- coreValues[[elem]]
            }
        }
    }
    slot(ggs, "params") <- params
    return(ggs)
}

#' Get number of functions from GenoGAMSetup
.nfun <- function(ggs) {
    vars <- .getVars(slot(ggs, "formula"), type = "covar")
    return(length(vars))
}

#' Get number of betas from GenoGAMSetup
.nbeta <- function(ggs) {
    betas <- slot(ggs, "beta")
    if(ncol(betas) > 1) {
        res <- nrow(betas)
    }
    else {
        funs <- .nfun(ggs)
        res <- nrow(betas)/funs
    }
    return(res)
}


#' Constructor function
#' @noRd
setupGenoGAM <- function(ggd, lambda = NULL, theta = NULL, H = 0, family = "nb", bpknots = 20, order = 2, penorder = 2) {

    ## knot placement    
    positions <- ranges(getIndex(ggd))[1]
    x <- start(positions):end(positions)
    nknots <- round(length(x)/bpknots)
    knots <- .placeKnots(x = x, nknots = nknots)

    X <- .buildDesignMatrix(ggd, knots = knots, pos = x, order = order)

    ## Number of betas = number of knots
    ## Number of functions = Count the functions in the formula
    nbetas <- nknots
    nfun <- length(.getVars(design(ggd), type = "covar"))
    S <- .buildSMatrix(nbetas, penorder, nfun)
    I <- .buildIMatrix(nbetas * nfun, H)
    S <- S + I

    ## turn knots into list to comply with object requirements
    knots <- list(knots)
    names(knots) <- "1"

    offset <- rep(sizeFactors(ggd), each = getTileSize(ggd))

    ggsetup <- GenoGAMSetup(params = list(lambda = lambda, theta = theta, H = H,
                                          order = order, penorder = penorder),
                            knots = knots, formula = design(ggd), 
                            offset = offset, family = family,
                            designMatrix = X, penaltyMatrix = S)
  
    return(ggsetup)
}

#' A function to generate knots for P-Splines from a GenoGAMDataSet object
#' @noRd
.generateKnotPositions <- function(ggd, bpknots = 20){

    positions <- IRanges::ranges(getIndex(ggd))[1]
    x <- IRanges::start(positions):IRanges::end(positions)
    nknots <- round(length(x)/bpknots)
    knots <- .placeKnots(x = x, nknots = nknots)

    res <- list(knots)
    names(res) <- "1"
    return(res)
}

#' A function to place knots for P-Splines
#' Courtesy to Simon Wood (mgcv). Slightly changed.
#' @noRd
.placeKnots <- function(x, nknots, ord = 2) {
  m <- ord + 1
  nk <- nknots - ord
  xu <- max(x)
  xl <- min(x)
  xr <- xu - xl
  multFactor <- min(1/(10^floor(log10(xr))), 0.001)
  xl <- xl - xr * multFactor
  xu <- xu + xr * multFactor
  dx <- (xu - xl)/(nk - 1)
  k <- seq(xl - dx * m, xu + dx * m, length.out = nk + 2 * ord + 2)
  
  return(k)
}

#' A function to build the penalization matrix S
#' Courtesy to Simon Wood (mgcv). Slightly changed
#' @noRd
.buildSMatrix <- function(p, order, nfun) {
    
    S <- Matrix::Matrix(diag(p), sparse = TRUE) ##initialize a diagonal identity matrix
    for (i in 1:order) {
        S <- diff(S) ## twice the difference   
    }
    S <- t(S)%*%S ## square
    design <- diag(nfun)
    res <- as(.blockMatrixFromDesignMatrix(S, design), "dgCMatrix")
    return(res)
}

#' A function to build the identity matrix I with multiple epsilon
#' Courtesy to Simon Wood (mgcv)
#' @noRd
.buildIMatrix <- function(p, epsilon) {
  I = Matrix::Matrix(diag(p), sparse = TRUE) ##initialize a diagonal identity matrix
  return(epsilon*I)
}

#B spline basis
.bspline <- function(x, k, ord = 2, derivative = 0) {
  res <- splines::spline.des(k, x, ord + 2, rep(derivative,length(x)), sparse=TRUE)$design
  return(res)
}

#' build a block matrix from a template submatrix and a design matrix
#' @noRd
.blockMatrixFromDesignMatrix <- function(template, design) {
  ## create 4-dim array by 'inserting' the template into the desing matrix
  arr <- array(template, c(dim(template),dim(design)))
  dims <- dim(arr)
  multP <- c(3,4,1,2)
  reduceP <- c(3,1,4,2)
  ## permute array for correct multiplication
  multArr <- aperm(arr, multP)*as.vector(design)
  ## permute array for correct reduction
  reducedArr <- aperm(multArr, reduceP)
  ## reduce 4-dim array to 2-dim matrix
  dim(reducedArr) <- c(nrow(template)*nrow(design), ncol(template)*ncol(design))
  return(reducedArr)
}

.buildDesignMatrix <- function(ggd, knots, pos, order) {

    ## build matrix
    x <- as(.bspline(pos, knots, order),"dgCMatrix")
    formulaCols <- .getVars(design(ggd))
    designCols <- as.vector(na.omit(formulaCols))
    design <- as.matrix(colData(ggd)[,designCols])
    if("s(x)" %in% names(formulaCols)) {
        control <- rep(1, nrow(design))
        design <- cbind(control, design)
    }
    X <- as(.blockMatrixFromDesignMatrix(x, design), "dgCMatrix")
    return(X)
}
