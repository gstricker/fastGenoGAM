% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readData.R
\name{readData}
\alias{readData}
\title{Read Data function}
\usage{
readData(config, hdf5 = FALSE, split = FALSE,
  settings = GenoGAMSettings(), ...)
}
\arguments{
\item{config}{A data.frame containing the experiment design of the model
to be computed with the first three columns fixed. See the 'experimentDesign'
parameter in \code{\link{GenoGAMDataSet}} or details here.}

\item{hdf5}{Should the data be stored on HDD in HDF5 format? By default this
is disabled, as the Rle representation of count data already provides a
decent compression of the data. However in case of large organisms, a complex
experiment design or just limited memory, this might further decrease the
memory footprint.}

\item{split}{If TRUE the data will be stored as a list of DataFrames by
chromosome instead of one big DataFrame. This is only necessary if organisms
with a genome size bigger than 2^31 (approx. 2.14Gbp) are analyzed,
in which case Rs lack of long integers prevents having a well compressed Rle
of sufficient size.}

\item{settings}{A GenoGAMSettings object. Not needed by default, but might
be of use if only specific regions should be read in.
See \code{\link{GenoGAMSettings}}.}

\item{...}{Further parameters that can be passed to low-level functions.
Mostly to pass arguments to custom process functions. In case the default
process functions are used, i.e. the default settings paramenter,
the most interesting parameters might be fragment length estimator method
from ?chipseq::estimate.mean.fraglen for single-end data.}
}
\value{
A DataFrame of counts for each sample and position.
Or if split = TRUE, a list of DataFrames by chromosomes
}
\description{
This is the core function to read and parse raw data from a config file.
At the moment only the BAM format is supported. It is not intended to be
used by the user directly, as it is called internally by the GenoGAMDataSet
constructor. However it is exported if people wish to separately assemble
their data and construct the GenoGAMDataSet from SummarizedExperiment
afterwards. It also offers the possibility to use the HDF5 backend.
}
\details{
The config data.frame contains the actual experiment design. It must
contain at least three columns with fixed names: 'ID', 'file' and 'paired'.

The field 'ID' stores a unique identifier for each alignment file.
It is recommended to use short and easy to understand identifiers because
they are subsequently used for labelling data and plots.

The field 'file' stores the complete path to the BAM file.

The field 'paired', values TRUE for paired-end sequencing data, and FALSE for
single-end sequencing data.

Other columns will be ignored by this function.
}
\examples{
# Read data

## Set config file
config <- system.file("extdata/Set1", "experimentDesign.txt", package = "fastGenoGAM")
config <- read.table(config, header = TRUE, sep = '\\t', stringsAsFactors = FALSE)
for(ii in 1:nrow(config)) {
    absPath <- system.file("extdata/Set1/bam", config$file[ii], package = "fastGenoGAM")
    config$file[ii] <- absPath
}

## Read all data
df <- readData(config)
df

## Read data of a particular chromosome
settings <- GenoGAMSettings(chromosomeList = "chrI")
df <- readData(config, settings = settings)
df

## Read data of particular range
region <- GenomicRanges::GRanges("chrI", IRanges(10000, 20000))
params <- Rsamtools::ScanBamParam(which = region)
settings <- GenoGAMSettings(bamParams = params)
df <- readData(config, settings = settings)
df
}
\author{
Georg Stricker \email{georg.stricker@in.tum.de}
}
