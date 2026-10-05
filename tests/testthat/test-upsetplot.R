library(epiSeeker)

context("test function for upsetplot")

test_that("upsetplot.csAnno works without vennpie", {
    skip_if_not_installed("ggupset")
    data("peakAnno", package = "epiSeeker")

    expect_s4_class(peakAnno, "csAnno")

    p <- upsetplot.csAnno(peakAnno, vennpie = FALSE)

    expect_true(
        inherits(p, "ggplot") ||
            inherits(p, "gtable") ||
            inherits(p, "patchwork")
    )
})


test_that("upsetplot.csAnno works with vennpie enabled", {
    skip_if_not_installed("ggupset")
    skip_if_not_installed("ggimage")
    data("peakAnno", package = "epiSeeker")

    expect_s4_class(peakAnno, "csAnno")

    p <- upsetplot.csAnno(peakAnno, vennpie = TRUE)

    expect_true(
        inherits(p, "ggplot") ||
            inherits(p, "gtable") ||
            inherits(p, "patchwork")
    )

    ## the plot has to be drawable, not only constructible: the vennpie
    ## sub-view is embedded as an annotation_custom() layer, which ggplot2 >= 4.0
    ## rejects outside of coord_cartesian()
    skip_if_not_installed("ggplot2")
    fn <- tempfile(fileext = ".png")
    on.exit(unlink(fn), add = TRUE)
    expect_error(
        suppressWarnings(ggplot2::ggsave(fn, p, width = 7, height = 5, dpi = 72)),
        NA
    )
})


test_that("upsetplot.csAnno with order_by", {
    skip_if_not_installed("ggupset")
    data("peakAnno", package = "epiSeeker")

    p <- upsetplot.csAnno(peakAnno, order_by = "degree")

    expect_true(
        inherits(p, "ggplot") ||
            inherits(p, "gtable") ||
            inherits(p, "patchwork")
    )
})
