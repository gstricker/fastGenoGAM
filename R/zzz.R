## clean up
.onUnload <- function (libpath) {
  library.dynam.unload("fastGenoGAM", libpath)
}
