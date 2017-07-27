## NOTE: All estimation and optimization is performed on the negative log-likelihood,
## as we are targeting a minimization problem (as is usually done in numerical optimization).
## The log-likelihood, the gradient and the hessian should usually be pre-multiplies by (-1)
## to get their normal form as derived on paper.

#' The IRLS algorithm
## .irls_nb <- function(beta0, X, y, theta, lambda, S, offset, control = list(eps = 1e-6, maxiter = 1000)) {
##     beta <- beta0
##     beta_old <- beta - control$eps - 1
##     dif <- max(abs(beta - beta_old))
##     k <- 0
##     converged <- TRUE
    
##     while (dif > control$eps & k < control$maxiter){
##         futile.logger::flog.debug(paste0("Iteration: ", k, "\n"))

##         eta <- offset + X %*% beta
##         mu <- exp(eta)
##         z <- 1/mu*(y - mu)+eta
  
##         ## 1. Compute V 
##         v <- as.numeric(mu^2/Matrix::diag(mu+(mu)^2/theta))
    
##         ## 2. Update weight matrix
##         W <- Matrix::bandSparse(dim(X)[1], k = 0, diag = c(list(v)))
        
##         ## 3. Compute new value of beta
##         beta_old <- beta
##         H <- Matrix::t(X) %*% W %*% X + 2*lambda * S
##         beta <- Matrix::solve(H, Matrix::t(X) %*% W %*% z)

##         dif <- max(abs(beta - beta_old))

##         futile.logger::flog.debug(paste0("Maximal Parameter Difference: ", dif, "\n"))
##         futile.logger::flog.debug(paste0("---------------------------\n"))
    
##         k <- k+1
##     }

##     if(k == control$maxiter) {
##         converged <- FALSE
##     }
##     res <- list(par = as.numeric(beta), converged = converged, iterations = k + 1)
##     return(res)
## }

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
        f <- .ll_penalized_nb
        gr <- .gr_ll_penalized_nb
        hf <- .hessian_nb
        args <- list(y = y, theta = params$theta)
    }

    
    H <- do.call(.compute_hessian, c(list(betas, X, offset, S, params$lambda, hf), args))

    ## fact argument needed to normalize likelihood with respect to data points
    res <- .lbfgs(x0 = betas, H0 = H, f = f, gr = gr, X = X, y = y,
                  offset = offset, theta = params$theta,
                  lambda = params$lambda, S = S, fact = dim(X)[1], control = control)

    res$par <- matrix(res$par, nrow = length(res$par), ncol = 1)
    return(res)
}

## NOTE: Likelihood and gradient have to have the same arguments, because
## arguments passed through the ellipse are passed to both functions. If
## they differ R will throw an error of unused parameter

#' penalized Negative Binomial likelihood
#' @noRd
.ll_penalized_nb <- function(beta, X, y, offset, theta, lambda, S, fact = 1){
  n <- dim(X)[1]
  eta <- offset + X%*%beta
  mu <- exp(eta)
  aux1 <- theta + y
  aux2 <- theta + mu
  ## pull the log inside to use gamma and factorial in log space due to 
  ## possibly very high numbers
  l <- sum(lgamma(aux1) - (lfactorial(y) + lgamma(theta))) + t(y) %*% eta + n*theta*log(theta) - t(aux1) %*% log(aux2)
  pen <- t(beta) %*% S %*% beta
  return(-1/fact * l[1] + lambda*pen[1,1])
}  

#' gradient of penalized negative Binomial likelihood
#' the gradients are divided by N to make lambda length independent
#' @noRd
.gr_ll_penalized_nb <- function(beta,X,y,offset,theta,lambda,S){
  eta <- offset + X%*%beta
  mu <- exp(eta)
  z <- (y-mu)/(1+mu/theta)
  res <- Matrix::t(X)%*%z
  pen <- S %*% beta
  return (-res[,1]+2*lambda*pen[,1])
}

#' Compute basepair Standard deviation from Hessian
#' @noRd
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
        f <- .hessian_nb
        args <- list(y = y, theta = params$theta)
    }
    H <- do.call(.compute_hessian, c(list(betas, X, offset, S, params$lambda, f), args))
    
    Hinv <- .invertHessian(H)

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
            se <- sqrt(Matrix::diag(Xsub%*%HinvSub%*%Matrix::t(Xsub)))
            return(se)
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
.compute_hessian <- function(beta, X, offset, S, lambda, f, ...){
    
    eta <- offset + X%*%beta    
    d <- f(eta, ...)
    
    if(length(d) == 0) {
        return(as(matrix(,0, 0), "dgCMatrix"))
    }
    
    D <- Matrix::bandSparse(dim(X)[1], k = 0, diag = c(list(d)))
    res <- (Matrix::t(X) %*% D %*% X + 2*lambda*S)
    return (res)
}

.hessian_nb <- function(eta, y, theta) {

    mu <- exp(eta)
    d <- mu * (1 + y/theta) / ((1 + mu/theta)^2)
    return(d)
}

## #' Computation of the inverse of H
## #' @noRd
## .invertHessian <- function(H, tol = 0.0001) {
##     if(length(H) == 0) {
##         return(as(matrix(, 0, 0), "dgCMatrix"))
##     }
##     Hinv <- invisible(Matrix::solve(H))
##     return(Hinv)
## }



#' Computation of the inverse of H
#' @noRd
.invertHessian <- function(H, tol = 0.0001) {

    if(length(H) == 0) {
        return(as(matrix(, 0, 0), "dgCMatrix"))
    }

    C <- Cholesky(H)
    
    res <- lapply(1:ncol(H), function(k) {
        ek = rep(0, nrow(H))
        ek[k] = 1
        s = Matrix::solve(C,ek)
        keep <- which(abs(s@x) > tol)
        x <- s@x[keep]
        i <- keep - 1
        j <- rep(k - 1, length(keep))
        return(list(x = x, i = i, j = j))
    })

    x <- unlist(sapply(res, function(m) m$x))
    i <- unlist(sapply(res, function(m) m$i))
    j <- unlist(sapply(res, function(m) m$j))

    if(length(x) == 0 | length(i) == 0 | length(j) == 0) {
        return(as(matrix(, 0, 0), "dgCMatrix"))
    }
    Hinv <- Matrix::sparseMatrix(i = i, j = j, x = x, index1 = FALSE)

    return(Hinv)
}

.computeDirection <- function(H, gradk, s, y, ro, len) {
    q <- gradk
    if(len == 0) {
        r <- Matrix::solve(H, q, sparse = TRUE)
        return(r)
    }

    alpha <- numeric(len)
    
    for(ii in len:1) {
        alpha[ii] <- ro[ii] * crossprod(s[,ii], q)
        q <- q - alpha[ii] * y[,ii]
    }

    r <- Matrix::solve(H, q, sparse = TRUE)
    
    for(jj in 1:len) {
        beta <- as.numeric(ro[jj] * crossprod(y[,jj], r))
        r <- r + s[,jj] * (alpha[jj] - beta)
    }
    return(r)
}

#' compute approximate Hessian
.Hupdate <- function(x, ...) {
    H <- .compute_hessian(beta = x, f = .hessian_nb, ...)
    return(H)
}

#' the simple backtrack algorithm for linesearch
#' to ensure Wolfe conditions. See Nocedal, Algorithm 3.1
.linesearch <- function(x, p, f, grad, alpha0 = 1, c = runif(1), rho = 0.9, ...) {
    alpha <- alpha0
    
    ## check if initial alpha is good enough
    betas <- as.numeric(x + alpha * p)
    fleft <- f(beta = betas, ...)
    fright <- as.numeric(f(beta = x, ...) + c * alpha * (Matrix::t(grad(beta = x, ...)) %*% p))

    ## decrease alpha by rho, to 'backtrack' it's optimal values
    while(fleft > fright) {
        alpha <- rho * alpha
        betas <- as.numeric(x + alpha * p)
        fleft <- f(beta = betas, ...)
        fright <- as.numeric(f(beta = x, ...) + c * alpha * (Matrix::t(grad(beta = x, ...)) %*% p))
    }
    return(alpha)
}

#' The modified lbfgs function for sparse optimization problems
#' @param x0 an initial vector
#' @param H0 an initial Hessian matrix. Can be missing, in which case
#' it will be the identity matrix
#' @param f the likelihood function
#' @param gr the likelihood gradient function
#' @param control parameters, if empty filled accordingly
#' @details Performs a modified L-BFGS. The modification stems largely from the fact,
#' that instead of an diagonal approximated inverse Hessian an approximated Hessian is
#' used and a linear system is solved via a direct solver is solved to obtain the
#' respective values. (see Nocedal, p. 178)
#' @noRd
.lbfgs <- function(x0, H0, f, gr, control = list(), fact = 1, ...) {
    ## If list is empty then replace with the proper settings
    if(length(control) == 0) {
        control <- list(eps = 1e-6, maxiter = 1000, alpha = 1, rho = 0.5, c = 1e-4, m = 6)
    }

    ## initialize Hessian as identity matrix
    if(missing(H0)) {
        H0 <- Matrix::bandSparse(length(x0), k = 0)
    }
    
    ## initialize variables
    k <- 1
    idx <- 0
    x <- x0
    H <- H0
    epsilon <- control$eps
    m <- control$m
    converged <- TRUE

    ## the only vectors stored.
    s <- matrix(nrow = length(x0), ncol = m)
    y <- matrix(nrow = length(x0), ncol = m)
    ro <- numeric(m)

    lllast <- 0
    llk <- epsilon + 1
    lldiff <- abs(llk - lllast)
    normGrad <- epsilon + 1
    grad <- matrix(gr(beta = x, ...), ncol = 1)
    
    ## Start L-BFGS loop
    while(lldiff > epsilon & normGrad > epsilon & k <= control$maxiter) {

        ## compute direction and stepsize
        p <- (-1) * .computeDirection(H, grad, s, y, ro, idx)
        
        alphak <- .linesearch(x, p, f, gr, alpha0 = control$alpha, c = control$c, rho = control$rho, ...)

        ## update params
        idx <- min(k, m) ## are we within the first m iterations?
        xnext <- matrix(x + alphak*p, ncol = 1)

        ## update ll and gradient
        lllast <- llk
        llk <- f(xnext, fact = fact, ...)
        lldiff <- abs(llk - lllast)
        gradNext <- matrix(gr(beta = xnext, ...), ncol = 1)
        normGrad <- sqrt(as.numeric(crossprod(gradNext)))

        ## update s, y and rho vectors
        ## if 
        if(k > m) {
            s <- cbind(s[,-1], 0)
            y <- cbind(y[,-1], 0)
            ro <- ro[-1]
        }

        s[,idx] <- xnext - x
        y[,idx] <- gradNext - grad
        ro[idx] <- as.numeric(1/crossprod(y[,idx], s[,idx]))

        

        ## create new Hessian
        H <- .Hupdate(xnext, ...)

        ## reset params variables
        x <- xnext
        grad <- gradNext

        ## printing
        futile.logger::flog.debug(paste0("Iteration: ", k, "\n"))
        futile.logger::flog.debug(paste0("||gr|| = ", normGrad, "\n"))
        futile.logger::flog.debug(paste0("Alpha = ", alphak, "\n"))
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




## setup <- fastGenoGAM:::.initiate(ggd, subggs, coords, id - 1)

## betas <- slot(setup, "beta")
## X <- slot(setup, "designMatrix")
## y <- slot(setup, "response")
## offset <- slot(setup, "offset")
## params <- slot(setup, "params")
## S <- slot(setup, "penaltyMatrix")

## library(profvis)

## slot(setup, "fits") <- fastGenoGAM:::.getFits(setup)
## H <- .compute_hessian_negbin(setup)

## futile.logger::flog.threshold(futile.logger::DEBUG)
## profvis(
##     res <- .lbfgs(x0 = betas, H0 = H, f = .ll_penalized_nb, gr = .gr_ll_penalized_nb, 
##                  X = X, y = y, offset = offset, theta = params$theta, lambda = params$lambda, S = S)
## )

