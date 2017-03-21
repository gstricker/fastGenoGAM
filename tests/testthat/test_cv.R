context("Testing Cross Validation functionality")

ggd <- makeTestGenoGAMDataSet(sim = TRUE)
ggs <- setupGenoGAM(ggd, lambda = NULL, theta = NULL, family = "nb",
                    H = 0, bpknots = 20, order = 2, penorder = 2)
folds <- 10
iv <- 20

test_that("the CV intervals are produced correctly", {
    cv <- .leaveOutConsecutiveIntervals(folds, iv, getTileSize(ggd)*2)
    expect_equal(length(cv), folds)
    expect_equal(length(cv[[1]]), getTileSize(ggd)/folds*2)
})

test_that("The likelihood function computes correctly", {
    id <- 3
    settings <- GenoGAMSettings()
    cv <- .leaveOutConsecutiveIntervals(folds, iv, getTileSize(ggd)*2)
    setup <- .initiate(ggd, ggs, id)
    fixedpars <- list(lambda = NULL, theta = NULL)
    initpars <- list(lambda = log(slot(setup, "params")[["lambda"]]),
                  theta = log(slot(setup, "params")[["theta"]]))

    ll <- .loglik(pars = initpars, setup = list(setup), CV_intervals = cv,
                  ov = getOverhangSize(ggd), method = slot(settings, "optimMethod"),
                  control = slot(settings, "optimControl"),
                  fixedpars = fixedpars)
    
    expect_true(ll < 0)
    expect_true(length(ll) == 1)

    fixedpars[] <- lapply(initpars, exp)
    ll_fixed <- .loglik(pars = initpars, setup = list(setup), CV_intervals = cv,
                  ov = getOverhangSize(ggd), method = slot(settings, "optimMethod"),
                  control = slot(settings, "optimControl"),
                  fixedpars = fixedpars)
    expect_true(all.equal(ll, ll_fixed))

    fixedpars$lambda <- NULL
    ll_fixed_one <- .loglik(pars = initpars, setup = list(setup), CV_intervals = cv,
                  ov = getOverhangSize(ggd), method = slot(settings, "optimMethod"),
                  control = slot(settings, "optimControl"),
                  fixedpars = fixedpars)
    expect_true(all.equal(ll, ll_fixed, ll_fixed_one))

    ll <- .loglik(pars = initpars, setup = list(setup), CV_intervals = cv,
                  ov = 0, method = slot(settings, "optimMethod"),
                  control = slot(settings, "optimControl"),
                  fixedpars = fixedpars)
    
    expect_true(ll < 0)
    expect_true(length(ll) == 1)

    ll <- .loglik(pars = initpars, setup = list(setup), CV_intervals = cv,
                  ov = getTileSize(ggd), method = slot(settings, "optimMethod"),
                  control = slot(settings, "optimControl"),
                  fixedpars = fixedpars)
    
    expect_true(ll == 0)
    expect_true(length(ll) == 1)
})
