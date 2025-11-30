library(epiSeeker)

context("test function for plotCov ")

test_that("plotCov works with sample peak", {
    
    peak <- readPeakFile(getSampleFiles()[[4]])

    p <- plotCov(peak = peak, weightCol = "V5")

    expect_true(inherits(p, "ggplot"))
})


test_that("getChrCov", {
    
    peak <- readPeakFile(getSampleFiles()[[4]])

    df <- getChrCov(
        peak = peak,
        weightCol = NULL,
        chrs = "chr1",
        xlim = c(0, 5e6)
    )

    expect_true(is.data.frame(df))
    expect_true(all(c("chr", "start", "end", "value") %in% names(df)))
    expect_true(nrow(df) >= 0)
})


test_that("sortChrName", {
    x <- c("chr10", "chr2", "chrX", "chr1")
    sorted <- sortChrName(x)

    expect_identical(sorted, c("chr1", "chr2", "chr10", "chrX"))
})
