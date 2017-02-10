context("Testing the GenoGAMDataSet-class")

## Unittests
gr <- GPos(GRanges("chr1", IRanges(1, 10000)))
seqlengths(gr) <- 1e6
df <- DataFrame(colA = 1:10000, colB = round(runif(10000)))
se <- SummarizedExperiment(rowRanges = gr, assays = list(df))
ggd <- GenoGAMDataSet(se, chunkSize = 2000, overhangSize = 250, 
                      design = ~ s(x) + s(x, by = "experiment"))

test_that("The constructor works correctly for SummarizedExperiment", {
    test_ggd <- ggd
    expect_identical(assay(test_ggd), assay(se))
    expect_identical(rowRanges(test_ggd), rowRanges(se))
    expect_is(test_ggd, "GenoGAMDataSet")

    gr <- GPos(GRanges("chr1", IRanges(1, 10000)))
    df <- DataFrame(colA = 1:10000, colB = round(runif(10000)))
    se <- SummarizedExperiment(rowRanges = gr, assays = list(df))
    expect_error(GenoGAMDataSet(se, chunkSize = 2000, overhangSize = 250, 
                                design = ~ s(x) + s(x, by = "experiment")),
                 "Sequence lengths missing in the Seqinfo object of SummarizedExperiment")

    gr <- GPos(GRanges(c("chr1", "chr1"), IRanges(c(1,5001), c(10000, 10000))))
    seqlengths(gr) <- 1e6
    df <- DataFrame(colA = 1:15000, colB = round(runif(15000)))
    se <- SummarizedExperiment(rowRanges = gr, assays = list(df))
    expect_error(GenoGAMDataSet(se, chunkSize = 2000, overhangSize = 250, 
                                design = ~ s(x) + s(x, by = "experiment")),
                 "Overlapping regions encountered. Please reduce ranges and data first.")
})

## test_that("The constructor works correctly for config files and data.frames", {
    
## })

test_that("Tiling works correctly under full chromosome conditions", {
    chromosomes <- GRanges(c("chr1", "chr2"), IRanges(c(1,1), c(10000, 10000)))
    seqlengths(chromosomes) <- c(1e6, 1e6)
    l <- list(chunkSize = 1000, overhangSize = 200, chromosomes = chromosomes)
    tiles <- .makeTiles(l)
    
    expect_equal(metadata(tiles)$numTiles, length(tiles))
    expect_equal(metadata(tiles)$numTiles, 20)
    expect_true(max(width(tiles)) == 1400 & min(width(tiles)) == 1400)
})

test_that("Tiling works correctly with no overhang", {
    chromosomes <- GRanges(c("chr1", "chr2"), IRanges(c(1,1), c(10000, 10000)))
    seqlengths(chromosomes) <- c(1e6, 1e6)
    l <- list(chunkSize = 1000, overhangSize = 0, chromosomes = chromosomes)
    tiles <- .makeTiles(l)
    
    expect_equal(metadata(tiles)$numTiles, length(tiles))
    expect_equal(metadata(tiles)$numTiles, 20)
    expect_true(max(width(tiles)) == 1000 & min(width(tiles)) == 1000)
})

test_that("Tiling works correctly with incorrect input", {
    chromosomes <- GRanges(c("chr1", "chr2"), IRanges(c(1,1), c(10000, 10000)))
    seqlengths(chromosomes) <- c(1e6, 1e6)
    l <- list(chunkSize = 1000, overhangSize = -100, chromosomes = chromosomes)
    expect_error(.makeTiles(l), "Overhang size must be equal or greater than 0")

    l <- list(chunkSize = -1000, overhangSize = 100, chromosomes = chromosomes)
    expect_error(.makeTiles(l), "Chunk size must be equal or greater than 1000")

    l <- list(chunkSize = 1000, overhangSize = 0, chromosomes = GRanges())
    expect_error(.makeTiles(l), "Chromosome list should contain at least one entry")
})

test_that("Tiling works correctly with a distinct set of regions", {
    chromosomes <- GRanges(c("chr1", "chr1", "chr2", "chr2"),
                           IRanges(c(1, 1200, 1, 10000), c(2000, 10000, 999, 11100)))
    seqlengths(chromosomes) <- c(1e6, 1e6)
    l <- list(chunkSize = 1000, overhangSize = 200, chromosomes = chromosomes)
    tiles <- .makeTiles(l)
    
    expect_equal(metadata(tiles)$numTiles, 12)
    expect_true(max(width(tiles)) == 1400 & min(width(tiles)) == 1400)
})

test_that("Settings checking functions work correctly", {
    test_ggd <- ggd

    expect_true(checkObject(test_ggd))
    
    metadata(slot(test_ggd, "index"))$check <- FALSE
    expect_false(checkObject(test_ggd))

    metadata(slot(test_ggd, "index"))$check <- NULL
    expect_false(checkObject(test_ggd))

    metadata(slot(test_ggd, "index"))$chunkSize <- 100
    expect_true(is(.checkChunkSize(test_ggd), "character"))

    metadata(slot(test_ggd, "index"))$tileSize <- 1000
    expect_true(is(.checkTileSize(test_ggd), "character"))

    end(slot(test_ggd, "index"))[1] <- 10000
    expect_true(is(.checkEqualityOfTiles(test_ggd), "character"))

    gr2 <- GRanges("chr2", IRanges(1, 1000), id = 30)
    seqlengths(gr2) <- 1e6
    md <- metadata(slot(test_ggd, "index"))
    suppressWarnings(slot(test_ggd, "index") <- c(slot(test_ggd, "index"), gr2))
    metadata(slot(test_ggd, "index")) <- md
    expect_true(is(.checkChromosomes(test_ggd), "character"))

    expect_true(is(.checkNumberOfTiles(test_ggd), "character"))

    metadata(slot(test_ggd, "index"))$chromosomes <- gr2
    expect_true(is(.checkTileRanges(test_ggd), "character"))
})

test_that("Accessors return the right slots", {
    test_ggd <- ggd

    expect_identical(getIndex(test_ggd), slot(test_ggd, "index"))
    expect_identical(tileSettings(test_ggd), metadata(slot(test_ggd, "index")))
    expect_identical(dataRange(test_ggd), rowRanges(test_ggd)@pos_runs)
    expect_identical(getChromosomes(test_ggd), metadata(slot(test_ggd, "index"))$chromosomes)
    expect_identical(getTileSize(test_ggd), metadata(slot(test_ggd, "index"))$tileSize)
    expect_identical(getChunkSize(test_ggd), metadata(slot(test_ggd, "index"))$chunkSize)
    expect_identical(getOverhangSize(test_ggd), metadata(slot(test_ggd, "index"))$overhangSize)
    expect_identical(getTileNumber(test_ggd), metadata(slot(test_ggd, "index"))$numTiles)
    expect_identical(design(test_ggd), slot(test_ggd, "design"))
    expect_identical(sizeFactors(test_ggd), slot(test_ggd, "sizeFactors"))
})


