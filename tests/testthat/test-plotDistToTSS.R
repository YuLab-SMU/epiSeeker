library(epiSeeker)

context("test function for plotDistToTSS ")

test_that("plotDistToTSS works with basic example", {
    data("peakAnno", package = "epiSeeker")

    df <- as.data.frame(peakAnno)


    p <- plotDistToTSS.data.frame(df, categoryColumn = 1)

    expect_s3_class(p, "ggplot")
})


test_that("plotDistToTSS works with custom distanceBreaks", {
    skip_if_not_installed("ggplot2")
    data("peakAnno", package = "epiSeeker")

    df <- as.data.frame(peakAnno)


    p <- plotDistToTSS.data.frame(
        df,
        distanceBreaks = c(0, 500, 2000, 5000, 20000),
        categoryColumn = 1
    )

    expect_s3_class(p, "ggplot")
})


test_that("plotDistToTSS works with palette option", {
    data("peakAnno", package = "epiSeeker")

    df <- as.data.frame(peakAnno)

    p <- plotDistToTSS.data.frame(
        df,
        palette = "Set3",
        categoryColumn = 1
    )

    expect_s3_class(p, "ggplot")
})


test_that("plotDistToTSS works with '.id'", {
    data("peakAnno", package = "epiSeeker")

    df <- as.data.frame(peakAnno)

    df$.id <- sample(c("A", "B"), size = nrow(df), replace = TRUE)

    expect_silent({
        p <- plotDistToTSS.data.frame(
            df,
            categoryColumn = ".id"
        )
    })

    expect_s3_class(p, "ggplot")
})
