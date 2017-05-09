library(devtools)
library(BiocCheck)
library(futile.logger)
flog.threshold(ERROR)
document("~/workspace/fastGenoGAM")
test("~/workspace/fastGenoGAM")
check("~/workspace/fastGenoGAM")
BiocCheck("~/workspace/fastGenoGAM")
install("~/workspace/fastGenoGAM")




######################
## Beta estimation tests
########################
library(fastGenoGAM)

## specify folder and experiment design path
folder <- system.file("extdata/Set1", package='GenoGAM')
expDesign <- file.path(folder, "experimentDesign.txt")

## specify chunk and overhang size, bpk being the knot spacing
bpknots <- 20
chunkSize <- 1000
overhangSize <- 15*bpknots

## build the GenoGAMDataSet
ggd <- GenoGAMDataSet(
  expDesign, directory = folder,
  chunkSize = chunkSize, overhangSize = overhangSize,
  design = ~ s(x) + s(x, by = genotype)
)

## compute size factors
ggd <- computeSizeFactors(ggd)

## restrict to the actual small data region (usually not needed, just for vignette purposes)
ggd <- subset(ggd, seqnames == "chrXIV" & pos >= 305001 & pos <= 308000)

lambda <- 4601
theta <- 4.51
family <- "nb"
H <- 0
order <- 2
m <- 2

## compute model without parameter estimation to save time in vignette
coords <- fastGenoGAM:::.getCoordinates(ggd)
chunks <- fastGenoGAM:::.getChunkCoords(coords)
ggs <- fastGenoGAM:::setupGenoGAM(ggd, lambda = lambda, theta = theta, family = family, 
                    H = H, bpknots = bpknots, order = order,
                    penorder = m)

id <- 2
control <- slot(slot(ggd, "settings"), "optimControl")
control$maxit <- control$betaMaxit
control$betaMaxit <- control$trace <- NULL
setup <- fastGenoGAM:::.initiate(ggd, ggs, coords, id)

##gold estimation
lambdaFun <- function(id, data, init, coords) {
    ## suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))

    control <- slot(slot(data, "settings"), "optimControl")
    control$maxit <- control$betaMaxit
    control$betaMaxit <- control$trace <- NULL
    setup <- fastGenoGAM:::.initiate(data, init, coords, id)
    betas <- fastGenoGAM:::.estimateParams(setup, control)
    print(betas$convergence)
    slot(setup, "beta") <- betas$par
    slot(setup, "fits") <- fastGenoGAM:::.getFits(setup)
     
    slot(setup, "designMatrix") <- new("dgCMatrix")
    slot(setup, "penaltyMatrix") <- new("dgCMatrix")
    slot(setup, "response") <- numeric()
    slot(setup, "offset") <- numeric()
    slot(setup, "params")$id <- id

    return(setup)
}

ids <- 1:length(coords)
res <- BiocParallel::bplapply(ids, lambdaFun, 
                                  data = ggd, init = ggs, coords = coords)
s <- start(chunks) %% start(coords) + 1
e <- s + width(chunks) - 1
relativeChunks <- IRanges::IRanges(s, e)

## combine fits
combinedFits <- fastGenoGAM:::.transformResults(res, relativeChunks, what = "fits")
png("converged.png", width = 1200, height = 800)
par(mfrow = c(3,1))
plot(assay(ggd)[,1], col = "#00000030")
lines(exp(combinedFits[,1]), col = "red")
plot(assay(ggd)[,3], col = "#00000030")
lines(exp(combinedFits[,2]), col = "red")
plot(assay(ggd)[,3], col = "#00000030")
lines(1:1600, exp(res[[1]]@fits[[2]][4801:6400]), col = "red")
lines(701:2300, exp(res[[2]]@fits[[2]][4801:6400]), col = "blue")
lines(1401:3000, exp(res[[3]]@fits[[2]][4801:6400]), col = "green")
dev.off()

## with maxcap
slot(slot(ggd, "settings"), "optimControl")$betaMaxit <- 100
resCapped <- BiocParallel::bplapply(ids, lambdaFun, 
                                  data = ggd, init = ggs, coords = coords)

## combine fits
combinedCappedFits <- fastGenoGAM:::.transformResults(resCapped, relativeChunks, what = "fits")
png("unconverged100.png", width = 1200, height = 800)
par(mfrow = c(3,1))
plot(assay(ggd)[,1], col = "#00000030")
lines(exp(combinedCappedFits[,1]), col = "red")
plot(assay(ggd)[,3], col = "#00000030")
lines(exp(combinedCappedFits[,2]), col = "red")
plot(assay(ggd)[,3], col = "#00000030")
lines(1:1600, exp(resCapped[[1]]@fits[[2]][4801:6400]), col = "red")
lines(701:2300, exp(resCapped[[2]]@fits[[2]][4801:6400]), col = "blue")
lines(1401:3000, exp(resCapped[[3]]@fits[[2]][4801:6400]), col = "green")
dev.off()

slot(slot(ggd, "settings"), "optimControl")$betaMaxit <- 2000
## estimation
goldTime <- system.time(gold <- fastGenoGAM:::.estimateParams(setup, control))

maxiter <- round(seq(100, gold$counts[1] - 100, length.out = 10))
allBetas <- NULL
allTimes <- NULL
allMSE <- NULL

for(ii in maxiter) {
    control$maxit <- ii
    time <- system.time(betas <- fastGenoGAM:::.estimateParams(setup, control))
    mse <- mean((betas$par - gold$par)^2)
    allMSE <- c(allMSE, mse)
    allBetas <- cbind(allBetas, betas$par)
    allTimes <- c(allTimes, time[["elapsed"]])
}

allBetas <- cbind(allBetas, gold$par)
allTimes <- c(allTimes, goldTime[["elapsed"]])
allMSE <- c(allMSE, 0)
iterations <- c(maxiter, gold$counts[1])

## plotting
plot(iterations, allTimes)
png("mse.png", width = 1200, height = 800)
plot(iterations, allMSE)
dev.off()

png("betas.png", width = 2200, height = 800)
par(mfrow = c(2,5))
for(ii in 1:length(maxiter)) {
    plot(allBetas[,ii], gold$par, xlab = maxiter[ii], cex.lab = 2, cex.axis = 2)
}
dev.off()

slot(setup, "beta") <- gold$par
fits <- as.data.frame(fastGenoGAM:::.getFits(setup))
fits <- fits[4801:nrow(fits),]
r <- coords[id,]
par(mfrow = c(2,1))
plot(assay(ggd)[start(r):end(r),1])
lines(exp(fits[,1]), col = "red")
plot(assay(ggd)[start(r):end(r),3])
lines(exp(fits[,2]), col = "red")
