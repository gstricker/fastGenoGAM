context("Testing the functionality of the main modelling functions")

ggd <- makeTestGenoGAMDataSet(sim = TRUE)
coords <- .getCoordinates(ggd)
emptyGGD <- GenoGAMDataSet()
emptyCoords <- .getCoordinates(emptyGGD)
ggs <- setupGenoGAM(ggd)
emptyGGS <- GenoGAMSetup()

test_that("The response vector is build correctly", {
    ## all combinations of different empty inputs
    res1 <- .buildResponseVector(emptyGGD, coords, 1)
    res2 <- .buildResponseVector(emptyGGD, emptyCoords, 1)
    res3 <- .buildResponseVector(ggd, emptyCoords, 1)
    
    expect_true(all(c(length(res1), length(res2), length(res3)) == 0))

    id <- sample(length(coords),1)
    res <- .buildResponseVector(ggd, coords, id)
    correctLength <- width(coords[id,]) * dim(ggd)[2]
    expect_true(length(res) == correctLength)
})

test_that("The specific tile setup initializes correctly", {
    ## all combinations of different empty inputs
    res1 <- .initiate(emptyGGD, ggs, coords, 1)
    params <- slot(res1, "params")
    expect_true(length(slot(res1, "response")) == 0)
    expect_true(all(is.na(slot(res1, "beta"))))
    expect_true(params$lambda == dim(slot(res1, "designMatrix"))[2])
    expect_true(params$theta == 1)

    res2 <- .initiate(ggd, emptyGGS, coords, 1)
    params <- slot(res2, "params")
    expect_true(length(slot(res1, "response")) == 0)
    expect_true(all(is.na(slot(res1, "beta"))))
    expect_true(params$lambda == 0)
    expect_true(params$theta == 0)

    res3 <- .initiate(ggd, emptyGGS, emptyCoords, 1)
    params <- slot(res3, "params")
    expect_true(length(slot(res1, "response")) == 0)
    expect_true(all(is.na(slot(res1, "beta"))))
    expect_true(params$lambda == 0)
    expect_true(params$theta == 0)

    res4 <- .initiate(emptyGGD, ggs, emptyCoords, 1)
    params <- slot(res1, "params")
    expect_true(length(slot(res1, "response")) == 0)
    expect_true(all(is.na(slot(res1, "beta"))))
    expect_true(params$lambda == dim(slot(res1, "designMatrix"))[2])
    expect_true(params$theta == 1)
    
    res5 <- .initiate(ggd, emptyGGS, emptyCoords, 1)
    params <- slot(res1, "params")
    expect_true(length(slot(res1, "response")) == 0)
    expect_true(all(is.na(slot(res1, "beta"))))
    expect_true(params$lambda == 0)
    expect_true(params$theta == 0)
})

test_that("Fits are correctly computed", {
    ## empty input
    res1 <- .getFits(emptyGGS)
    expect_true(length(res1) == 0)

    ## only empty betas
    res2 <- .initiate(emptyGGD, ggs, coords, 1)
    expect_true(length(res1) == 0)

    ## non-empty input with betas == 1, for easier checking
    ## fits should be just the row sums of the design matrix
    ## which should be 1.
    test_ggs <- ggs
    slot(test_ggs, "beta") <- matrix(1, dim(slot(test_ggs, "designMatrix"))[2], 1)
    res3 <- .getFits(test_ggs)
    X <- as.matrix(slot(test_ggs, "designMatrix"))
    fits1 <- rowSums(X[1:(dim(X)[1]/2), 1:(dim(X)[2]/2)])
    fits2 <- rowSums(X[(dim(X)[1]/2 + 1):dim(X)[1], (dim(X)[2]/2 + 1):dim(X)[2]])
    default <- rep(1, dim(X)[2]/2)
    expect_true(length(res3) == 2)
    expect_true(all.equal(res3[[1]], fits1, default))
    expect_true(all.equal(res3[[2]], fits2, default))
    fitNames <- c("s(x)", paste("s(x)", colnames(colData(ggd)), sep = ":"))
    expect_true(all(names(res3) == fitNames))
})

## NEXT TEST ESTIMATEPARAMS
