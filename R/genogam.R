#################################
## The main modelling function ##
#################################

##' genogam
##'
##' The main modelling function.
##' @param ggd The GenoGAMDataSet object to be fitted
##' @param lambda The penalization parameter. If NULL (default) estimated by
##' cross validation. So far only one parameter for all splines is supported.
##' @param theta The global overdispersion parameter. If NULL (default) estimated
##' by cross validation.
##' @param family The distribution to be used. So far only Negative-Binomial (nb)
##' is supported.
##' @param H The factor for additional first-order regularization. This should be
##' zero (default) in most cases. It can be useful for sparse data with many
##' zero-counts regions or very low coverage. In this cases it is advised to use
##' a small factor like 0.01, which would penalize those regions but not the ones
##' with higher coverage. See Wood S., Generalized Additive Models (2006) for more.
##' @param bpknots The spacing between knots. The lower the number, the more local
##' functions and the more sensitive the model with the trade-off of increased
##' computation time. In our experience going below 20 is not necessary even for
##' close double peaks.
##' @param kfolds The number of folds for cross validation
##' @param intervalSize The size of the hold-out intervals in cross validation.
##' If replicates are present, this can easily be increased to twice the fragment
##' size to capture more of the local correlation. If no replicates are present,
##' keep the number low to avoid heavy interpolation (default).
##' @param regions How many regions should be used in cross validation? The number
##' is an upper limit. It is usually corrected down, such that the total number of
##' models computed during cross validation does not exceed the total number of
##' models to compute for the entire genome. This is usually the case for small
##' organisms such as yeast.
##' @param order The order of the B-spline basis, which is order + 2, where 0
##' is the lowest order. Thus order = 2 is equivalent to cubic order (= 3).
##' @param m The order of penalization. Thus m = 2 penalizes the second differences.
##' @return The fit as a GenoGAM object
##' @examples
##' ggd <- makeTestGenoGAMDataSet(sim = TRUE)
##' res <- genogam(ggd, lambda = 266.8368, theta = 2.415738)
##' @author Georg Stricker \email{georg.stricker@@in.tum.de}
##' @export
genogam <- function(ggd, lambda = NULL, theta = NULL, family = "nb", H = 0,
                    bpknots = 20, kfolds = 10, intervalSize = 20, 
                    regions = 20, regionSize = 4000, order = 2, m = 2) {

    futile.logger::flog.info("Initializing the model")

    input <- paste0("Initializing the model with the following parameters:\n",
                    "  Lambda: ", lambda, "\n",
                    "  Theta: ", theta, "\n",
                    "  Family: ", family, "\n",
                    "  H: ", H, "\n",
                    "  Knot spacing: ", bpknots, "\n",
                    "  Number of folds: ", kfolds, "\n",
                    "  Interval size: ", intervalSize, "\n",
                    "  Number of regions: ", regions, "\n",
                    "  B-spline order: ", order, "\n",
                    "  Penalization order: ", m, "\n",
                    "  GenoGAMDataSet object: ")
    futile.logger::flog.debug(input)
    futile.logger::flog.debug(show(ggd))
    futile.logger::flog.debug(show(slot(ggd, "settings")))
    
    settings <- slot(ggd, "settings")
    check <- .checkSettings(ggd)
    coords <- .getCoordinates(ggd)
    chunks <- .getChunkCoords(coords)

    if(!check) break

    ggs <- setupGenoGAM(ggd, lambda = lambda, theta = theta, family = family, 
                        H = H, bpknots = bpknots, order = order,
                        penorder = m, control = slot(settings, "estimControl"))

    futile.logger::flog.info("Done")
    ## Cross Validation
    cv <- FALSE

    if(is.null(lambda) | is.null(theta)) {
        futile.logger::flog.info("Estimating hyperparameters")

        ## check correct region size
        if(regionSize < 2*getOverhangSize(ggd)) {
            futile.logger::flog.warning("region size too small for cross validation, set to twice the overhang size.")
            regionSize <- 2*getOverhangSize(ggd)
        }

        ## backup original tile index
        index_backup <- getIndex(ggd)

        ## make new tile index
        metadata(slot(ggd, "index"))$chunkSize <- regionSize
        newTiles <- .makeTiles(tileSettings(ggd))
        slot(ggd, "index") <- newTiles
        new_ggs <- setupGenoGAM(ggd, lambda = lambda, theta = theta, family = family, 
                                H = H, bpknots = bpknots, order = order,
                                penorder = m, control = slot(settings, "estimControl"))
        new_coords <- .getCoordinates(ggd)
        
        ## get the tile ids for CV
        sumMatrix <- sum(ggd)
        ## ncv set to a number, such that CV does not compute more models than
        ## the actual genogam run. 
        ## 40 is the average expected number of iterations
        ncv <- min(regions, ceiling(length(new_coords)/(kfolds*40)))

        if(ncv < regions) {
            futile.logger::flog.debug(paste("Reducing the number of regions to", ncv))
        }
        
        if(ncv < length(new_coords)) {
            if(sum(sapply(colData(ggd), sum)) == nrow(colData(ggd)) |
               nrow(sumMatrix) < regions) {
                rsums <- rowSums(sumMatrix)
                ids <- order(rsums, decreasing = TRUE)[1:ncv]
            }
            else {
                futile.logger::flog.debug("Performing DESeq2 analysis to find regions with highest fold change")
                pvals <- suppressMessages(suppressWarnings(.deseq(sumMatrix, colData(ggd))))
                ids <- order(pvals)[1:ncv]
            }
        }
        else ids <- 1:length(new_coords)

        control <- slot(settings, "optimControl")
        futile.logger::flog.debug(paste("Selected the following regions:", paste(ids, collapse = ",")))
        
        params <- .doCrossValidation(ggd, setup = new_ggs, coords = new_coords, 
                                     id = ids, folds = kfolds, 
                                     intervalSize = intervalSize,
                                     fn = .loglik, 
                                     method = slot(settings, "optimMethod"),
                                     control = control) 

        names(params) <- c("lambda", "theta")
        slot(ggs, "params")$lambda <- params["lambda"]
        slot(ggs, "params")$theta <- params["theta"]

        ## set back the original tile index and delete other objects
        slot(ggd, "index") <- index_backup
        rm("new_ggs", "new_coords")
        gc()
        
        cv <- TRUE
        futile.logger::flog.info("Done")
    }
    
    futile.logger::flog.info("Fitting model")

    .local <- function(id, data, init, coords) {
        ## suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))

        setup <- .initiate(data, init, coords, id)
        betas <- .estimateParams(setup)

        futile.logger::flog.debug(paste("Beta estimation for tile", id, "took", betas$iterations, "iterations"))
        if(betas$converged == FALSE) {
            futile.logger::flog.warn("Beta estimation did not converge. Increasing the 'maxiter' or 'eps' parameter in 'estimControl' slot in the settings might help, but should be done at own risk.")
        }
        
        slot(setup, "beta") <- betas$par
        slot(setup, "fits") <- .getFits(setup)
        
        slot(setup, "se") <- .compute_SE(setup)
        
        slot(setup, "designMatrix") <- new("dgCMatrix")
        slot(setup, "penaltyMatrix") <- new("dgCMatrix")
        slot(setup, "response") <- numeric()
        slot(setup, "offset") <- numeric()
        slot(setup, "params")$id <- id

        return(setup)
    }

    ids <- 1:length(coords)
    res <- BiocParallel::bplapply(ids, .local, 
                                  data = ggd, init = ggs, coords = coords)

    futile.logger::flog.info("Done")

    futile.logger::flog.info("Assembling fits and building GenoGAM object")
    ## make chunk coordinates relative
    s <- start(chunks) %% start(coords) + 1
    e <- s + width(chunks) - 1
    relativeChunks <- IRanges::IRanges(s, e)

    ## combine fits
    combinedFits <- .transformResults(res, relativeChunks, what = "fits")
    combinedSEs <- .transformResults(res, relativeChunks, what = "se")
        
    ## build GenoGAM object
    se <- SummarizedExperiment::SummarizedExperiment(rowRanges = SummarizedExperiment::rowRanges(ggd),
                                                     assays = list(fits = combinedFits,
                                                                   se = combinedSEs))
    
    modelParams <- slot(ggs, "params")
    modelParams$cv <- cv
    modelParams$bpknots <- bpknots
    modelParams$order <- order
    modelParams$penorder <- m
    if(cv) {
        modelParams$kfolds <- kfolds
        modelParams$intervalSize <- intervalSize
        modelParams$regions <- ncv
    }
    modelParams <- c(modelParams, tileSettings(ggd))
    modelParams$check <- NULL

    gg <- GenoGAM(se, family = slot(ggs, "family"),
                   design = design(ggd),
                   sizeFactors = sizeFactors(ggd),
                   factorialDesign = colData(ggd),
                   params = modelParams,
                   settings = settings) ## smooths slot to be added
    
    ## How to make it run and write directly parallel. Do parallel by Chromosome, write to HDD and run next chromosome?

    futile.logger::flog.info("Done")
    return(gg)
}

#' Builds the response vector from GenoGAMDataSet and the given row coordinates
#' @noRd
.buildResponseVector <- function(ggd, coords, id) {
    
    if(length(ggd) == 0 | length(coords) == 0) {
        return(numeric())
    }
    
    tile <- coords[id,]
    y <- .subsetByCoords(x = assay(ggd), i = IRanges::start(tile):IRanges::end(tile))
    
    Y <- unname(unlist(as.data.frame(y)))
    return(as.integer(Y))
}

#' initiates GenoGAMSetup with tile specific data
#' @noRd
.initiate <- function(ggd, setup, coords, id) {

    ## initiate response vector and betas
    slot(setup, "response") <- .buildResponseVector(ggd, coords, id)
    numBetas <- dim(slot(setup, "designMatrix"))[2]
    
    if(length(slot(setup, "response")) != 0) {
        ## set betas to the runmedian of the response
        means <- IRanges::runmed(slot(setup, "response"), 11, endrule = "median")
        ## take 250 values at equidistant positions
        ks <- as.integer(seq(1, length(means), length.out = numBetas))
        ## set all zero values to 1, because of log
        betas <- means[ks]
        betas[which(betas == 0)] <- 1
        slot(setup, "beta") <- matrix(log(betas), numBetas, 1)
    }
        
    ## initialize lambda and theta
    params <- slot(setup, "params")
    if(is.null(params$lambda)) {
        params$lambda <- numBetas
    }
    if(is.null(params$theta)) {
        params$theta <- 1
    }
    slot(setup, "params") <- params

    return(setup)
}

#' compute fits from design matrix and estimated betas
#' @noRd
.getFits <- function(setup, log = TRUE) {
    if(all(dim(slot(setup, "designMatrix")) == c(0, 0)) |
       all(is.na(slot(setup, "beta")))) {
        return(list())
    }

    ## get design matrix and beta vector
    X <- slot(setup, "designMatrix")
    betas <- slot(setup, "beta")

    ## get row and column indeces of the design matrix
    ## for the respective fits
    des <- slot(setup, "design")
    nSplines <- ncol(des)
    nSamples <- nrow(des)
    dims <- dim(X)
    rows <- list(1:dims[1])
    cols <- list(1:dims[2])

    if(nSplines > 1) {
        cols <- split(1:dims[2], cut(1:dims[2], nSplines, labels = 1:nSplines))
    }
    if(nSamples > 1) {
        rows <- split(1:dims[1], cut(1:dims[1], nSamples, labels = 1:nSamples))
    }

    ## compute the fits by column (splines) and row (samples)
    ## of the design matrix. All fits of the same spline
    ## are combined to a vector to have the same structure as
    ## the response. Each spline is one element of the list, that
    ## gets returned as the result.
    fits <- lapply(1:nSplines, function(j) {
        res <- sapply(1:nSamples, function(i) {
            Xsub <- X[rows[[i]], cols[[j]]]
            beta <- betas[cols[[j]], 1]
            fit <- as.vector(Xsub %*% beta)
            if(!log) {
                fit <- exp(fit)
            }
            return(fit)
        })
        res <- as.vector(res)
        return(res)
    })
    varNames <- .makeNames(slot(setup, "formula"))
    names(fits) <- varNames

    return(fits)
}

##' transform the result list of GenoGAMSetup object into a DataFrame
##' of either fits or standard errors
##' @noRd
.transformResults <- function(x, relativeChunks, what = c("fits", "se")) {
    if(length(relativeChunks) == 0) {
        return(data.table::data.table())
    }

    allData <- lapply(x, function(y) {
        if(length(slot(y, what)) == 0) {
            return(data.table::data.table())
        }
            
        s <- start(relativeChunks[slot(y, "params")$id])
        e <- end(relativeChunks[slot(y, "params")$id])

        ## go by column and subset by colData column
        l <- lapply(1:length(slot(y, what)), function(z) {
            
            des <- slot(y, "design")
            fits <- slot(y, "fits")
            if(length(des) == 0 | length(fits) == 0) {
                res <- list()
            }
            else {
                tileLength <- length(fits[[z]])/nrow(des)
                ## expand colData columns to full length and coerce to boolean
                desVec <- matrix(as.logical(rep(des, each = tileLength)), ncol = ncol(des))
                ## take the z list element (= column) and
                ## subset by the z colData column
                res <- slot(y, what)[[z]][desVec[,z]]
                res <- res[s:e]
            }
            return(res)
        })
        names(l) <- names(slot(y, what))
        return(data.table::as.data.table(l))
    })

    combinedData <- data.table::rbindlist(allData)
    varNames <- colnames(combinedData)
    combinedData <- S4Vectors::DataFrame(combinedData)
    colnames(combinedData) <- varNames

    return(combinedData)
}
