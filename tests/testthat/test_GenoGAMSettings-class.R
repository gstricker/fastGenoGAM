context("Test the GenoGAMSettings-class constructor")

## Unittests
test_that("The GenoGAMSettings constructor works correctly", {
    ggs <- new("GenoGAMSettings")
    test_ggs <- GenoGAMSettings()
    expect_identical(ggs, test_ggs)
})

test_that("The GenoGAMSettings validation works correctly", {
    expect_error(GenoGAMSettings(bamParams = "myParams"))
    expect_error(GenoGAMSettings(optimMethod = 10))
    expect_error(GenoGAMSettings(optimControl = 10))
    expect_error(GenoGAMSettings(center = 10))
    expect_error(GenoGAMSettings(chromosomeList = 10))
    expect_error(GenoGAMSettings(processFunction = 10))
})

test_that("The GenoGAMSettings assignment works correctly", {
    params <- Rsamtools::ScanBamParam(what = c("pos"))
    ggs <- GenoGAMSettings(bamParams = params)
    expect_identical(params, slot(ggs, "bamParams"))
    
    params <- "myOptimMethod"
    ggs <- GenoGAMSettings(optimMethod = params)
    expect_identical(params, slot(ggs, "optimMethod"))

    params <- list(a = 1, b = "2")
    ggs <- GenoGAMSettings(optimControl = params)
    expect_identical(params, slot(ggs, "optimControl"))
    
    params <- TRUE
    ggs <- GenoGAMSettings(center = params)
    expect_identical(params, slot(ggs, "center"))

    params <- c("chr1", "chr2")
    ggs <- GenoGAMSettings(chromosomeList = params)
    expect_identical(params, slot(ggs, "chromosomeList"))

    params <- identity
    ggs <- GenoGAMSettings(processFunction = params)
    expect_identical(params, slot(ggs, "processFunction"))
})
