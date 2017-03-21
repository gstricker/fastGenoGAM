#################################
## The main modelling function ##
#################################

## genogam <- function(ggd, lambda = NULL, theta = NULL, family = "nb", H = 0,
##                     bpknots = 20, kfolds = 10, intervalSize = 20, 
##                     intervals = 20, order = 2, m = 2) {
    
##     settings <- slot(ggd, "settings")
##     tileSettings <- tileSettings(ggd)
##     check <- .checkSettings(ggd)
##     coords <- .getCoordinates(ggd)

##     if(!check) break

##     ggs <- setupGenoGAM(ggd, lambda = lambda, theta = theta, family = family, 
##                         H = H, bpknots = bpknots, order = order,
##                         penorder = m)

##     ## Cross Validation
##     cv <- FALSE

##     if(is.null(lambda) | is.null(theta)) {
##         futile.logger::flog.info("Estimating parameters")
        
##         ## get the tile ids for CV
##         sumMatrix <- sum(ggd)
##         ## ncv set to a number, such that CV does not compute more models than
##         ## the actual genogam run. 
##         ## 50 is the average expected number of iterations
##         ncv <- min(intervals, ceiling(length(coords)/(kfolds*50)))
##         if(ncv < length(coords)) {
##             if(sum(sapply(colData(ggd), sum)) == nrow(colData(ggd))) {
##                 rsums <- rowSums(sumMatrix[[1]])
##                 ids <- order(rsums, decreasing = TRUE)[1:ncv]
##             }
##             else {
##                 pvals <- suppressMessages(suppressWarnings(.deseq(sumMatrix, colData(ggd))))
##                 ids <- order(pvals)[1:ncv]
##             }
##         }
##         else ids <- 1:length(coords)
    
##         params <- .doCrossValidation(ggd, setup = ggs, coords = coords, 
##                                      id = ids, folds = kfolds, 
##                                      intervalSize = intervalSize,
##                                      fn = .loglik, 
##                                      method = "Nelder-Mead", #GenoGAM:::getDefaults(settings, "optimMethod"),
##                                      control = GenoGAM:::getDefaults(settings, "optimControl")) 
##         proc.time() - start
##         names(params) <- c("lambda", "theta")
##         lambda <- params[1]
##         theta <- params[2]
##         slot(ggs, "params")$lambda <- lambda
##         slot(ggs, "params")$theta <- theta
        
##         cv <- TRUE
##         futile.logger::flog.info("Done")
##     }
    
##     futile.logger::flog.info("Fitting model")

##     lambdaFun <- function(id, data, init, coords) {
##         ## suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))

##         setup <- .initiate(data, init, coords, id)
##         betas <- .estimateParams(setup)
##         slot(setup, "beta") <- betas$par
##         slot(setup, "fits") <- .getFits(setup)
##         slot(setup, "designMatrix") <- new("dgCMatrix")
##         slot(setup, "penaltyMatrix") <- new("dgCMatrix")
##         slot(setup, "response") <- numeric()

##         return(setup)
##     }

##     ids <- coords$id
##     res <- BiocParallel::bplapply(ids, lambdaFun, 
##                                   data = ggd, init = ggs, coords = coords)
## }

.buildResponseVector <- function(ggd, id) {
    index <- getIndex(ggd)
    tile <- index[index$id == id,]
    positions <- rowRanges(ggd)
    rows <- subjectHits(findOverlaps(tile, positions))
    y <- assay(ggd)[rows,]
    
    Y <- unname(unlist(as.data.frame(y)))
    return(Y)
}


.initiate <- function(ggd, setup, id) {

    ## initiate response vector and betas
    slot(setup, "response") <- .buildResponseVector(ggd, id)
    numBetas <- dim(slot(setup, "designMatrix"))[2]
    slot(setup, "beta") <- matrix(mean(slot(setup, "response"), na.rm = TRUE), numBetas, 1)
        
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
                 method = "L-BFGS", control = list(fnscale=-1, maxit = 1000))
    return(res)
}

#penalized likelihood to be maximized \ negbin
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

#gradient of penalized likelihood \negbin
.gradient_likelihood_penalized_nb <- function(beta,X,y,offset,theta,lambda,S){
  eta <- offset + X%*%beta
  mu <- exp(eta)
  z <- (y-mu)/(1+mu/theta)
  res <- t(X)%*%z
  pen <- S %*% beta
  return (res[,1]-2*lambda*pen[,1])
}
