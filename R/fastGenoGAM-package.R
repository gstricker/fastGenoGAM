#' fastGenoGAM: A package providing a framework to analyse ChIP-Seq data
#' 
# @name fastGenoGAM
#' @import methods
#' @import HDF5Array
#' @import rhdf5
#' @import BiocParallel
#' @import IRanges
#' @import GenomicRanges
#' @importFrom stats runif rnbinom as.formula dnbinom optim na.omit
#' @importFrom futile.logger flog.info
#' @importFrom futile.logger flog.warn
#' @importFrom futile.logger flog.error
#' @importFrom futile.logger flog.trace
#' @importFrom futile.logger flog.debug
#' @importFrom futile.logger flog.threshold
#' @useDynLib fastGenoGAM
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
