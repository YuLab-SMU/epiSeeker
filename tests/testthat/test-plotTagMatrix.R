library(epiSeeker)

context("test function for plot tagMatrix")

test_that("plotPeakHeatmap_sub works with single tagMatrix", {
    data("tagMatrix", package = "epiSeeker")

    p <- plotPeakHeatmap_sub(tagMatrix = tagMatrix)


    expect_true(inherits(p, "ggplot"))
})



test_that("plotPeakHeatmap_sub works with list tagMatrix", {
    data("tagMatrix", package = "epiSeeker")

    tm_list <- list(a = tagMatrix, b = tagMatrix)

    p <- plotPeakHeatmap_sub(tm_list)

    expect_true(
        inherits(p, "ggplot") ||
        inherits(p, "gtable") ||
        inherits(p, "patchwork")
    )
})



test_that("plotPeakProf works with tagMatrix", {
    data("tagMatrix", package = "epiSeeker")

    p <- plotPeakProf(tagMatrix)

    expect_true(
        inherits(p, "ggplot")
    )
})



test_that("plotPeakHeatmap works to plot heatmap and profile)", {
    data("tagMatrix", package = "epiSeeker")

    p <- plotPeakHeatmap(tagMatrix, plot_prof = TRUE)

    expect_true(
        inherits(p, "ggplot") ||
        inherits(p, "gtable") ||
        inherits(p, "patchwork")
    )
})



test_that("plotPeakHeatmap works without profile", {
    data("tagMatrix", package = "epiSeeker")

    p <- plotPeakHeatmap(tagMatrix, plot_prof = FALSE)

    expect_true(
        inherits(p, "ggplot") ||
        inherits(p, "gtable") ||
        inherits(p, "patchwork")
    )
})
