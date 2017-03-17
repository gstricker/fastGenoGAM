## context("Testing Cross Validation functionality")

## ggd <- makeTestGenoGAMDataSet(sim = TRUE)
## ggs <- setupGenoGAM(ggd, lambda = NULL, theta = NULL, family = "nb",
##                     H = 0, bpknots = 20, order = 2, penorder = 2)
## test_that("the CV intervals are produced correctly", {
##     folds <- 10
##     iv <- 20
##     cv <- .leaveOutConsecutiveIntervals(folds, iv, 2000)
##     expect_equal(length(cv), folds)
##     expect_equal(length(cv[[1]]), folds * iv)
## })

## test_that("The likelihood function computes correctly", {
    
## })

