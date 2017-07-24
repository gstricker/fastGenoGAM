#' get varibles from a formula object
#' @noRd
.getVars <- function(formula, type = c("by", "covar")) {
    type <- match.arg(type)
    chFormula <- as.character(formula)
    variables <- gsub(" ", "", strsplit(chFormula[length(chFormula)],"\\+")[[1]])
    r <- switch(type,
                by = "by=(.+)\\)",
                covar = "s\\((.+)\\)")
    if(type == "covar") variables <- gsub(",.+", ")", variables)
    temp <- sapply(sapply(variables, function(m) regmatches(m,regexec(r,m))), function(y) y[2])
    res <- gsub(",(.+)", "", temp)
    return(res)
}

.makeNames <- function(formula) {
    vars <- .getVars(formula)
    vars[!is.na(vars)] <- paste("s(x)", vars[!is.na(vars)], sep = ":")
    vars[is.na(vars)] <- "s(x)"
    return(vars)
}

#' helper function to help fill a list of parameters
#' with the defaults if some are missing
.fillParameters <- function(l, ...) {
    params <- c(...)
    allin <- names(params) %in% names(l)

    ## has to be done after the 'allin', since if all params in l
    ## are not valid, %in% would match an empty list and not fill
    ## it with the default values
    wrong_params <- names(l) %in% names(params)
    if(sum(!wrong_params) > 0) {
        futile.logger::flog.warn("Some supplied parameters aren't valid and won't be used")
        l <- l[wrong_params]
    }
    
    if(!all(allin)) {
        for(elem in names(params)) {
            if(is.null(l[[elem]])) {
                l[[elem]] <- params[[elem]]
            }
        }
    }
    return(l)
}

## #' #' Update Formula with a specific penalization parameter lambda
## #' Not used at the moment
## #'
## #' @param formula A formula object
## #' @param lambda A vector of lambda values
## #' @return Updated formula object
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## .updateFormula <- function(formula, knots, m) {
##     covars <- .getVars(formula, "covar")
##     chFormula <- as.character(formula)
##     variables <- strsplit(chFormula[2],"\\+")[[1]]
##     variables[covars == "x"] <- gsub("x(?=,|\\))", "pos", variables[covars == "x"], perl = TRUE)
       
##     variables <- gsub(")", paste0(', bs = "ps", k = ', knots ,', m = ', m ,')'),
##                                 variables)
##     variablesCombined <- paste0(variables,collapse = " + ")
##     reassembledFormula <- as.formula(paste("value ~ offset(offset) + ", variablesCombined))
##     return(reassembledFormula)
## }

## #' Subset the object according to certain position in the coords.
## #' 
## #' @param object A data.frame alike object.
## #' @param A vector of length two, specifying the start and end sites to be trimmed to
## #' @return A trimmed object
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## .untile <- function(object, coords = NULL){
##     if(is.null(coords)) coords <- c(start(metadata(object)$chunks), end(metadata(object)$chunks))
##     start <- min(object$pos)
##     end <- max(object$pos)
##     overhang <- c(which(object$pos < coords[1]), which(object$pos >= coords[2]))
##     return(object[-overhang,])
## }

## #' Create a melted version of GenomicTiles,adding the experimentMatrix columns
## #' and melting by samplenames.
## #'
## #' @param gtiles, A DataFrameList alike object.
## #' @param settings A list of tile settings.
## #' @param colData The sample specific values, to be added to the melted data.frame.
## #' The rownames of colData should comply with the sample names in the config data.frame.
## #' @param experimentMatrix A matrix object representing the experiment design.
## #' @return A melted data.frame with additional columns.
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## .meltGTile <- function(gtiles, experimentDesign, sf, formula) {
##     res <- lapply(gtiles, function(gt) {

##         ids <- rownames(experimentDesign)
##         ## melt
##         melted <- reshape2::melt(as.data.frame(gt), measure.vars = ids,
##              variable.name = "ID")

##         ## add offset
##         if(!("offset" %in% names(melted))) {
##             melted$offset <- 0
##         }

##         rle <- Rle(melted$ID)
##         sf <- sf[match(runValue(rle), names(sf))] ## put them in right order
##         melted$offset <- melted$offset + rep(sf, runLength(rle))
        
##         ## add experimentDesign columns
##         designVars <- .getVars(formula, "by")
##         designCols <- experimentDesign[names(sf),na.omit(designVars), drop = FALSE]
##         for(col in names(designCols)) {
##             melted[[col]] <- rep(designCols[[col]], runLength(rle))
##         }
##         return(melted)
##     })
##     names(res) <- names(gtiles)
##     return(res)
## }
    
