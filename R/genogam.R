#################################
## The main modelling function ##
#################################

## lambda <- 4608.962 
## theta <- 1.504998 
## ## lambda = NULL
## ## theta = NULL
## family = "nb"
## H = 0
## bpknots = 20
## kfolds = 10
## intervalSize = 20
## intervals = 20
## order = 2
## m = 2

genogam <- function(ggd, lambda = NULL, theta = NULL, family = "nb", H = 0,
                    bpknots = 20, kfolds = 10, intervalSize = 20, 
                    intervals = 20, order = 2, m = 2) {
    
    settings <- slot(ggd, "settings")
    tileSettings <- tileSettings(ggd)
    check <- .checkSettings(ggd)
    coords <- .getCoordinates(ggd)
    chunks <- .getChunkCoords(coords)

    if(!check) break

    ggs <- setupGenoGAM(ggd, lambda = lambda, theta = theta, family = family, 
                        H = H, bpknots = bpknots, order = order,
                        penorder = m)

    ## Cross Validation
    cv <- FALSE

    if(is.null(lambda) | is.null(theta)) {
        futile.logger::flog.info("Estimating parameters")
        
        ## get the tile ids for CV
        sumMatrix <- sum(ggd)
        ## ncv set to a number, such that CV does not compute more models than
        ## the actual genogam run. 
        ## 50 is the average expected number of iterations
        ncv <- min(intervals, ceiling(length(coords)/(kfolds*50)))
        if(ncv < length(coords)) {
            if(sum(sapply(colData(ggd), sum)) == nrow(colData(ggd))) {
                rsums <- rowSums(sumMatrix[[1]])
                ids <- order(rsums, decreasing = TRUE)[1:ncv]
            }
            else {
                pvals <- suppressMessages(suppressWarnings(.deseq(sumMatrix, colData(ggd))))
                ids <- order(pvals)[1:ncv]
            }
        }
        else ids <- 1:length(coords)
    
        params <- .doCrossValidation(ggd, setup = ggs, coords = coords, 
                                     id = ids, folds = kfolds, 
                                     intervalSize = intervalSize,
                                     fn = .loglik, 
                                     method = slot(settings, "optimMethod"),
                                     control = slot(settings, "optimControl")) 

        names(params) <- c("lambda", "theta")
        slot(ggs, "params")$lambda <- params["lambda"]
        slot(ggs, "params")$theta <- params["theta"]
                
        cv <- TRUE
        futile.logger::flog.info("Done")
    }
    
    futile.logger::flog.info("Fitting model")

    lambdaFun <- function(id, data, init, coords) {
        ## suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))

        setup <- .initiate(data, init, coords, id)
        betas <- .estimateParams(setup)
        slot(setup, "beta") <- betas$par
        slot(setup, "fits") <- .getFits(setup)
        
        H <- .compute_hessian_negbin(setup)
        slot(setup, "sd") <- .computeSD(setup, H)
        
        slot(setup, "designMatrix") <- new("dgCMatrix")
        slot(setup, "penaltyMatrix") <- new("dgCMatrix")
        slot(setup, "response") <- numeric()
        slot(setup, "offset") <- numeric()
        slot(setup, "params")$id <- id

        return(setup)
    }

    ids <- 1:length(coords)
    res <- BiocParallel::bplapply(ids, lambdaFun, 
                                  data = ggd, init = ggs, coords = coords)

    ## make chunk coordinates relative
    s <- start(chunks) %% start(coords) + 1
    e <- s + width(chunks) - 1
    relativeChunks <- IRanges::IRanges(s, e)

    ## combine fits
    allFits <- lapply(res, function(y) {
        s <- start(relativeChunks[slot(y, "params")$id])
        e <- end(relativeChunks[slot(y, "params")$id])
        dt <- cbind(data.table::as.data.table(slot(y, "fits"))[s:e,], data.table::as.data.table(slot(y, "sd"))[s:e,])
        return(dt)
    })

    combined <- data.table::rbindlist(allFits)
    varNames <- colnames(combined)
    combined <- DataFrame(combined)
    colnames(combined) <- varNames

    ## build GenoGAM object
    se <- SummarizedExperiment::SummarizedExperiment(rowRanges = rowRanges(ggd),
                                                     assays = list(combined))
    modelParams <- slot(ggs, "params")
    modelParams$cv <- cv
    modelParams$bpknots <- bpknots
    modelParams$order <- order
    modelParams$penorder <- m
    if(cv) {
        modelParams$kfolds <- kfolds
        modelParams$intervalSize <- intervalSize
        modelParams$intervals <- ncv
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

    return(gg)
}

#' Builds the response vector from GenoGAMDataSet and the given row coordinates
#' @noRd
.buildResponseVector <- function(ggd, coords, id) {
    tile <- coords[id,]
    y <- assay(ggd)[IRanges::start(tile):IRanges::end(tile),]
    
    Y <- unname(unlist(as.data.frame(y)))
    return(Y)
}

#' initiates GenoGAMSetup with tile specific data
#' @noRd
.initiate <- function(ggd, setup, coords, id) {

    ## initiate response vector and betas
    slot(setup, "response") <- .buildResponseVector(ggd, coords, id)
    numBetas <- dim(slot(setup, "designMatrix"))[2]
    slot(setup, "beta") <- matrix(log(mean(slot(setup, "response"), na.rm = TRUE)), numBetas, 1)
        
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
    X <- slot(setup, "designMatrix")
    betas <- slot(setup, "beta")
    
    nSplines <- .nfun(setup)
    dims <- dim(X)
    if(nSplines > 1) {
        rows <- split(1:dims[1], cut(1:dims[1], nSplines, labels = 1:nSplines))
        cols <- split(1:dims[2], cut(1:dims[2], nSplines, labels = 1:nSplines))
    }
    else {
        rows <- list(1:dims[1])
        cols <- list(1:dims[2])
    }
    fits <- lapply(1:nSplines, function(y) {
        Xsub <- X[rows[[y]], cols[[y]]]
        beta <- betas[cols[[y]], 1]
        res <- as.vector(Xsub %*% beta)
        if(!log) {
            res <- exp(res)
        }
        return(res)
    })
    varNames <- .makeNames(slot(setup, "formula"))
    names(fits) <- varNames

    return(fits)
}

#' estimate betas
#' @noRd
.estimateParams <- function(ggs) {
    ## turn the beta matrix into a 1-column matrix
    betas <- slot(ggs, "beta")

    distr <- slot(ggs, "family")
    X <- slot(ggs, "designMatrix")
    y <- slot(ggs, "response")
    offset <- slot(ggs, "offset")
    params <- slot(ggs, "params")
    S <- slot(ggs, "penaltyMatrix")

    if (distr == "nb") {    
        likelihood <- .likelihood_penalized_nb
        gradient <- .gradient_likelihood_penalized_nb
    }

    res <- optim(betas, likelihood, gradient, X = X, y = y, offset = offset,
                 theta = params$theta, lambda = params$lambda, S = S, 
                 method = "L-BFGS", control = list(fnscale=-1, maxit = 200))
    return(res)
}

#' penalized Negative Binomial likelihood
#' @noRd
.likelihood_penalized_nb <- function(beta,X,y,offset,theta,lambda,S){
  n <- dim(X)[1]
  eta <- offset + X%*%beta
  mu <- exp(eta)
  aux1 <- theta + y
  aux2 <- theta + mu
  ## pull the log inside to use gamma and factorial in log space due to 
  ## possibly very high numbers
  l <- sum(lgamma(aux1) - (lfactorial(y) + lgamma(theta))) + t(y) %*% eta + n*theta*log(theta) - t(aux1) %*% log(aux2)
  pen <- t(beta) %*% S %*% beta
  return(l[1]-lambda*pen[1,1])
}  

#' gradient of penalized negative Binomial likelihood
#' @noRd
.gradient_likelihood_penalized_nb <- function(beta,X,y,offset,theta,lambda,S){
  eta <- offset + X%*%beta
  mu <- exp(eta)
  z <- (y-mu)/(1+mu/theta)
  res <- t(X)%*%z
  pen <- S %*% beta
  return (res[,1]-2*lambda*pen[,1])
}

#' Compute the penalized Hessian for Negative Binomial
#' @noRd
.compute_hessian_negbin <- function(setup){
    params <- slot(setup, "params")
    theta <- params$theta
    lambda <- params$lambda
    offset <- slot(setup, "offset")
    fits <- slot(setup, "fits")
    X <- slot(setup, "designMatrix")
    y <- slot(setup, "response")
    S <- slot(setup, "penaltyMatrix")
    
    eta <- offset + rowSums(as.data.frame(fits))
    names(eta) <- NULL
    
    mu <- exp(eta)
    d <- - mu * (1 + y/theta) / ((1 + mu/theta)^2) 
    res <- t(X) %*% Matrix::bandSparse(dim(X)[1], k = 0, diag = c(list(d)))%*% X - 2*lambda*S 
    return (res)
}

#' Computation of the inverse of H
#' @noRd
.invertHessian <- function(H, tol = 0.0001) {
    Hinv <- Matrix::solve(H)
    Hinv[abs(Hinv) < tol] <- 0
    return(Hinv)
}

#' Compute basepair Standard deviation from Hessian
#' @noRd
.computeSD <- function(setup, H) {
    X <- slot(setup, "designMatrix")
    Hinv <- .invertHessian(H)*(-1)
    
    nSplines <- .nfun(setup)
    dims <- dim(X)
    if(nSplines > 1) {
        rows <- split(1:dims[1], cut(1:dims[1], nSplines, labels = 1:nSplines))
        cols <- split(1:dims[2], cut(1:dims[2], nSplines, labels = 1:nSplines))
    }
    else {
        rows <- list(1:dims[1])
        cols <- list(1:dims[2])
    }
    sds <- lapply(1:nSplines, function(y) {
        Xsub <- X[rows[[y]], cols[[y]]]
        HinvSub <- Hinv[cols[[y]],cols[[y]]]
        res <- sqrt(Matrix::diag(Xsub%*%HinvSub%*%t(Xsub)))
        return(res)
    })
    varNames <- .makeNames(slot(setup, "formula"))
    names(sds) <- paste("se", varNames, sep = ".")

    return(sds)
}
