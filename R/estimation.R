## NOTE: All estimation and optimization is performed on the negative log-likelihood,
## as we are targeting a minimization problem (as is usually done in numerical optimization).
## The log-likelihood, the gradient and the hessian should usually be pre-multiplies by (-1)
## to get their normal form as derived on paper.

#' estimate betas
#' @noRd
.estimateParams <- function(ggs) {

    betas <- slot(ggs, "beta")

    distr <- slot(ggs, "family")
    X <- slot(ggs, "designMatrix")
    y <- slot(ggs, "response")
    offset <- slot(ggs, "offset")
    params <- slot(ggs, "params")
    S <- slot(ggs, "penaltyMatrix")
    control <- slot(ggs, "control")

    if(length(ggs) == 0) {
        res <- list(par = matrix(0,0,0), converged = FALSE, iterations = 0)
        return(res)
    }

    if (distr == "nb") {
        f <- ll_pen_nb
        gr <- gr_ll_pen_nb
        args <- list(beta = betas, X = X, XT = Matrix::t(X), offset = unname(offset),
                     y = y, S = S, lambda = params$lambda, theta = params$theta)
    }

    if(sum(slot(ggs, "response")) == 0){
        par <- rep(0, nrow(betas))
        res <- list(par = matrix(par, nrow = length(par), ncol = 1), converged = TRUE, iterations = 0)
        matrix(res$par, nrow = length(res$par), ncol = 1)
    }
    else {
        H <- do.call(.compute_hessian, c(list(family = distr), args))

        res <- .newton(x0 = betas, H0 = H, f = f, gr = gr, X = X, y = y,
                       offset = offset, theta = params$theta,
                       lambda = params$lambda, S = S, fact = dim(X), control = control)
    }
        
    return(res)
}

.compute_SE <- function(setup){
    if(length(setup) == 0) {
        return(numeric())
    }
    
    params <- slot(setup, "params")
    theta <- params$theta
    lambda <- params$lambda
    offset <- slot(setup, "offset")
    betas <- slot(setup, "beta")
    X <- slot(setup, "designMatrix")
    y <- slot(setup, "response")
    S <- slot(setup, "penaltyMatrix")
    distr <- slot(setup, "family")
    des <- slot(setup, "design")

    if(distr == "nb") {
        args <- list(beta = betas, X = X, XT = Matrix::t(X), offset = unname(offset),
                     y = y, S = S, lambda = params$lambda, theta = params$theta)
    }

    H <- do.call(.compute_hessian, c(list(family = distr), args))

    batchsize <- slot(setup, "control")$batchsize
    Hinv <- .invertHessian(H, batchsize = batchsize)

    ## get indeces of the design matrix and the Hessian
    ## by row and column according to the design
    
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

    ## compute the standard error by column (splines) and row (samples)
    ## of the design matrix and the Hessian. All standard errors of the
    ## same spline are combined to a vector to have the same structure as
    ## the response. Each spline is one element of the list, that
    ## gets returned as the result.
    ses <- lapply(1:nSplines, function(j) {
        res <- sapply(1:nSamples, function(i) {
            Xsub <- X[rows[[i]], cols[[j]]]
            HinvSub <- Hinv[cols[[j]],cols[[j]]]
            se <- compute_stdError(Xsub, HinvSub)
            if(length(se@x) != 0) {
                return(se@x)
            }
            else {
                return(rep(0L, nrow(Xsub)))
            }
        })
        res <- as.vector(res)
        return(res)
    })
    varNames <- .makeNames(slot(setup, "formula"))
    names(ses) <- varNames

    return(ses)
}

#' Compute the penalized Hessian
#' @noRd
.compute_hessian <- function(family, ...){
    res <- as(matrix(,0, 0), "dgCMatrix")
    
    if(family == "nb") {
        res <- compute_pen_hessian(...)
    }
    
    ## if(length(d) == 0) {
    ##     return(as(matrix(,0, 0), "dgCMatrix"))
    ## }
    
    return(res)
}

#' Computation of the inverse of H
#' @noRd
.invertHessian <- function(H, batchsize = 100) {
    if(length(H) == 0) {
        return(as(matrix(, 0, 0), "dgCMatrix"))
    }
        
    ## finding optimal batch size, is taken care of in the setup
    ## if it fails here look at the setupGenoGAM function
    nbatches <- ncol(H)/batchsize
    
    ij <- summary(H)
    ij <- ij[,-3]
    ij$jj <- (ij - 1)$j%/%batchsize + 1
    ij$ii <- ij$i + dim(H)[1]*(ij$j - 1) - batchsize*dim(H)[1]*(ij$jj - 1)
    pointers <- split(ij$ii, ij$jj)

    C <- Matrix::Cholesky(H)
    
    res <- lapply(1:nbatches, function(k) {
        top <- matrix(0, batchsize*(k - 1), batchsize)
        middle <- diag(1, batchsize, batchsize)
        bottom <- matrix(0, dim(H)[1] - batchsize*k, batchsize)
        ek <- rbind(top, middle, bottom)
        s = Matrix::solve(C,ek)
        keep <- pointers[[k]]
        x <- s@x[keep]
        return(x)
    })

    x <- unlist(res)

    Hinv <- Matrix::sparseMatrix(i = ij$i, j = ij$j, x = x, index1 = TRUE)

    return(Hinv)
}

.Hupdate <- function(x, X = X, XT = XT, ...) {
    
    H <- compute_pen_hessian(beta = x, X = X, XT = XT, ...)
    return(H)
    
}

.newton <- function(x0, H0, f, gr, X, control = list(), fact = c(1,1), ...) {
    ## If list is empty then replace with the proper settings
    if(length(control) == 0) {
        control <- list(eps = 1e-6, maxiter = 1000)
    }

    ## initialize Hessian as identity matrix
    if(missing(H0)) {
        H0 <- Matrix::bandSparse(length(x0), k = 0)
    }
    
    ## initialize variables
    k <- 1
    x <- x0
    H <- H0
    XT <- Matrix::t(X)
    epsilon <- control$eps
    converged <- TRUE

    lllast <- 0
    llk <- epsilon + 1
    lldiff <- abs(llk - lllast)
    normGrad <- epsilon + 1
    grad <- gr(beta = x, X = X, XT = XT, ...)
    
    ## Start L-BFGS loop
    while(lldiff > epsilon & normGrad > epsilon & k <= control$maxiter) {

        ## compute direction
        p <- (-1) * Matrix::solve(H, grad)@x
        
        ## update params
        xnext <- matrix(x + p, ncol = 1)

        ## update ll and gradient
        lllast <- llk
        llk <- f(xnext, X = X, ll_factor = fact[1], lambda_factor = fact[2], n = fact[1], ...)
        lldiff <- abs(llk - lllast)
        gradNext <- gr(beta = xnext, X = X, XT = XT, ...)
        normGrad <- sqrt(as.numeric(crossprod(gradNext))) 

        ## create new Hessian
        H <- .Hupdate(xnext, X = X, XT = XT, ...)

        ## reset params variables
        x <- xnext
        grad <- gradNext

        ## printing
        futile.logger::flog.debug(paste0("Iteration: ", k, "\n"))
        futile.logger::flog.debug(paste0("||gr|| = ", normGrad, "\n"))
        futile.logger::flog.debug(paste0("Likelihood = ", llk, "\n"))
        futile.logger::flog.debug("---------------------------\n")

        ## update iteration
        k <- k + 1
    }

    if(k == control$maxiter) {
        converged <- FALSE
    }
    return(list(par = x, converged = converged, iterations = k))
}

#' finds all possible batch sizes based on prime number factorization
.findBatchsize <- function(x) {
    facs <- factorize(x)
    res <- unlist(sapply(1:length(facs), function(y) {
        combinations <- combn(facs, y)
        apply(combinations, 2, prod)
    }))
    return(res)
}

#' Copied unchanged from Bill Venables package 'conf.design'
#' Version: 2.0.0
#' Published: 2013-02-23 15:18:29
#' 
#' Slightly changed the output
factorize <- function(x, divisors = primes(max(x)), ...) {
  stopifnot(is.numeric(x))
  if (length(x) > 1) {
    l <- vector("list", length(x))
    names(l) <- as.character(x)
    for (i in seq_along(x))
      l[[i]] <- Recall(x[i], divisors = divisors, ...)
    return(l)
  }
  if (x %% 1 > 0 || x < 2)
    return(x)
  tab <- divisors
  fac <- numeric(0)
  while(length(tab <- tab[x %% tab == 0]) > 0) {
    x <- x/prod(tab)
    fac <- c(fac, tab)
  }
  
  sort(fac)
}

#' Copied unchanged from Bill Venables package 'conf.design'
#' Version: 2.0.0
#' Published: 2013-02-23 15:18:29
primes <- function(n) {
### Find all primes less than n (or max(n) if length(n) > 1).
### Uses an obvious sieve method.  Nothing flash.
###
### 2013: This now uses a slightly improved coding of the version of
###       the algorithm used in the pracma package primes() function
###
  stopifnot(is.numeric(n), all(n %% 1 == 0))
  if ((M2 <- max(n)) <= 1)
    return(numeric(0))
  x <- seq(1, M2, 2)
  np <- length(x)
  x[1] <- 2
  if(M2 > 8) {
    top <- floor(sqrt(M2))
    p <- 1
    while((p <- p+2) <= top)
        if(x[(p + 1)/2] > 0)
            x[seq((p*p + 1)/2, np, p)] <- 0
  }
  x[x > 0]
}
