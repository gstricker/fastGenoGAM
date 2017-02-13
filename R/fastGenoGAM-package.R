#' fastGenoGAM: A package providing a framework to analyse ChIP-Seq data
#' 
#' @name fastGenoGAM
#' @import S4Vectors
#' @import BiocParallel
#' @import IRanges
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @import HDF5Array
#' @import methods
#' @importFrom Rsamtools bamWhich
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamHeader
#' @importFrom data.table data.table
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#' @importFrom data.table setnames
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomeInfoDb seqlevelsInUse
#' @importFrom futile.logger flog.info
#' @importFrom futile.logger flog.warn
#' @importFrom futile.logger flog.error
#' @importFrom futile.logger flog.trace
#' @importFrom futile.logger flog.debug
#' @import DESeq2
#' @import Biostrings
"_PACKAGE"
