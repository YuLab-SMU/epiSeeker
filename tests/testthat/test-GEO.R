library(epiSeeker)

context("test function for GEO data mining")


test_that("get_gsminfo", {
    gs <- get_gsminfo()
    expect_true(is.data.frame(gs))
})
