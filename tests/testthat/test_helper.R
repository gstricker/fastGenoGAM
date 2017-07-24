context("Test the correctness of some helper functions")

## Unittests
test_that("The parameter validation is correct", {
    l <- list(a = 1, b = 2, d = "a")
    params <- list(a = 5, b = 4, c = 10)
    true_params <- list(a = 1, b = 2, c = 10)
    validParams <- .fillParameters(l, params)
    valid <- sapply(1:length(validParams), function(ii) {
        true_params[[ii]] == validParams[[ii]]
    })
    expect_true(all(valid))
})
