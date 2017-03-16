#################################
## The main modelling function ##
#################################

genogam <- function(ggd, lambda = NULL, theta = NULL, family = "nb", H = 0,
                    bpknots = 20, kfolds = 10, intervalSize = 20, 
                    intervals = 20, order = 2, m = 2) {
    
    settings <- slot(ggd, "settings")
    tileSettings <- tileSettings(ggd)
    check <- .checkSettings(ggd)
    coords <- .getCoordinates(ggd)

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
                                     fn = .loglik, ov = getOverhangSize(ggd), 
                                     method = "Nelder-Mead", #GenoGAM:::getDefaults(settings, "optimMethod"),
                                     control = GenoGAM:::getDefaults(settings, "optimControl")) 
        proc.time() - start
        names(params) <- c("lambda", "theta")
        lambda <- params[1]
        theta <- params[2]
        slot(ggs, "params")$lambda <- lambda
        slot(ggs, "params")$theta <- theta
        
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
        slot(setup, "designMatrix") <- new("dgCMatrix")
        slot(setup, "penaltyMatrix") <- new("dgCMatrix")
        slot(setup, "response") <- numeric()

        return(setup)
    }

    ids <- coords$id
    res <- BiocParallel::bplapply(ids, lambdaFun, 
                                  data = ggd, init = ggs, coords = coords)
}

