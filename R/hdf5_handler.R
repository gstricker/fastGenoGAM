##' make HDF5 filename for initial pre-processed data
.makeFilename <- function(dir, name) {
    date <- strsplit(as.character(Sys.time()), " ")[[1]][1]
    suffix <- gsub("-", "", date)
    f <- file.path(dir, paste0(name, "_", suffix, ".h5"))
    return(f)
}

##' Function to write DataFrame pre-processed counts to HDF5
##' @noRd
.writeToHDF5 <- function(df, file, name = file, settings, simple = FALSE) {
  
    if(futile.logger::flog.threshold() == "DEBUG") {
        verbose <- TRUE
    }
    else {
        verbose <- FALSE
    }

    dir <- slot(settings, "hdf5Control")$dir
    if(!dir.exists(dir)) {
        futile.logger::flog.info(paste("HDF5 directory created at:", dir))
        dir.create(dir)
    }

    if(!is.null(slot(settings, "hdf5Control")$chunk)) {
        chunkdims <- slot(settings, "hdf5Control")$chunk
    }
    else {
        chunkdims <- HDF5Array::getHDF5DumpChunkDim(dim(df), "integer")
    }

    futile.logger::flog.info(paste("Writing", name, "to HDF5"))
    if(simple) {
        h5file <- file.path(dir, file)
    }
    else {
        h5file <- .makeFilename(dir, file)
    }
    if(file.exists(h5file)) {
        stop(paste("File", h5file, "exists. Please remove before a new file can be written."))
    }
    h5 <- HDF5Array::writeHDF5Array(HDF5Array::DelayedArray(df), file = h5file, name = name, chunk_dim = chunkdims, verbose = verbose)
    futile.logger::flog.info(paste(name, "written"))
    return(h5)
}


##' create HDF5 file
.createH5File <- function(seed, dir, what = c("fits", "coefs")) {
    what <- match.arg(what)

    if(what == "coefs") {
        ident <- strsplit(seed, split = "_")[[1]]
    }
    if(what == "fits") {
        ident <- strsplit(seed, split = "/")[[1]]
    }
    
    suffix <- ident[length(ident)]
    f <- paste0(what, "_", suffix)
    h5file <- rhdf5::H5Fcreate(file.path(dir, f))
    return(list(pointer = h5file, file = file.path(dir, f)))
}

##' initialize HDF5 dataset with dimensions d
.createH5DF <- function(ggd, d, id = 1, what = c("fits", "coefs")) {   
    what <- match.arg(what)

    seedFile <- assay(ggd)[[id]]@seed@file
    dir <- ggd@settings@hdf5Control$dir

    h5file <- .createH5File(seedFile, dir, what)
   
    h5space <- rhdf5::H5Screate_simple(d,d)
    if(what == "fits") {
        h5fits <- rhdf5::H5Dcreate(h5file$pointer, "fits", "H5T_IEEE_F32LE", h5space) ## Type Float: "H5T_IEEE_F32LE"
        h5ses <- rhdf5::H5Dcreate(h5file$pointer, "ses", "H5T_IEEE_F32LE", h5space)
        H5Dclose(h5ses)
    }
    if(what == "coefs") {
        h5fits <- rhdf5::H5Dcreate(h5file$pointer, "coefs", "H5T_IEEE_F32LE", h5space)
    }

    H5Dclose(h5fits)
    H5Sclose(h5space)
    H5Fclose(h5file$pointer)
    return(h5file$file)
}

### The folowing is modified from HDF5Array package. Thanks to Herve Pagès.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A simple mechanism to lock/unlock a file so processes can get temporary
### exclusive access to it
###

##' creates the lock file
.locked_path <- function(filepath)
{
    if (!isSingleString(filepath) || filepath == "") {
        stop("'filepath' must be a single non-empty string")
    }
    dir <- dirname(filepath)
    return(file.path(dir, ".lock"))
}

##' locks the file by creating an empty ".lock" file
.lock_file <- function(filepath)
{
    locked_path <- .locked_path(filepath)
    
    ## Must wait if the file is already locked.
    while (file.exists(locked_path)) {
        Sys.sleep(runif(1, 0, 0.5))
    }

    if(!file.create(locked_path)) {
        stop("failed to lock '", filepath, "' file")
    }
    return(locked_path)
}

##' unlocks file by removing the empty ".lock" file
.unlock_file <- function(lock)
{
    if (!file.remove(lock)) {
        stop("failed to unlock '", filepath, "' file")
    }
}
