library(epiSeeker)

context("test function for upsetplot")

test_that("upsetplot.csAnno works without vennpie", {
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
    data("peakAnno", package = "epiSeeker")

    expect_s4_class(peakAnno, "csAnno")

    p <- upsetplot.csAnno(peakAnno, vennpie = TRUE)

    expect_true(
        inherits(p, "ggplot") ||
            inherits(p, "gtable") ||
            inherits(p, "patchwork")
    )
})


test_that("upsetplot.csAnno with order_by", {
    data("peakAnno", package = "epiSeeker")

    p <- upsetplot.csAnno(peakAnno, order_by = "degree")

    expect_true(
        inherits(p, "ggplot") ||
            inherits(p, "gtable") ||
            inherits(p, "patchwork")
    )
})
