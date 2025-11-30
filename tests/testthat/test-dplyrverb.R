library(epiSeeker)

context("test function for dplyr verb ")

test_that("filter.GRanges for demo_peak", {

    data("demo_peak", package = "epiSeeker")

    df <- as.data.frame(demo_peak)
    col_for_filter <- names(df)[ncol(df)]

    res <- filter(demo_peak, !!as.name(col_for_filter) > median(df[[col_for_filter]]))

    expect_s4_class(res, "GRanges")
    expect_true(length(res) <= length(demo_peak))
})


test_that("mutate.GRanges for demo_peak", {

    data("demo_peak", package = "epiSeeker")

    res1 <- dplyr::mutate(demo_peak, new_var = 1)
    expect_s4_class(res1, "GRanges")
    expect_true("new_var" %in% names(mcols(res1)))

    res2 <- dplyr::mutate(demo_peak, new2 = 2, .before = 1)
    expect_s4_class(res2, "GRanges")

    res3 <- dplyr::mutate(demo_peak, new3 = 3, .after = 1)
    expect_s4_class(res3, "GRanges")
})


test_that("rename.GRanges", {

    data("demo_peak", package = "epiSeeker")

    df <- as.data.frame(demo_peak)
    old_col <- names(df)[ncol(df)]  
    new_col <- "renamed_col"

    res <- rename(demo_peak, !!new_col := !!as.name(old_col))

    expect_s4_class(res, "GRanges")
    expect_true(new_col %in% names(mcols(res)))
})


test_that("arrange.GRanges", {

    data("demo_peak", package = "epiSeeker")

    res <- dplyr::arrange(demo_peak, seqnames)

    expect_s4_class(res, "GRanges")
    expect_equal(length(res), length(demo_peak))
})
