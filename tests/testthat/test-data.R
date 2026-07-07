library(epiSeeker)

context("test build-in data")

test_that("demo_peak", {
    data("demo_peak", package = "epiSeeker")

    expect_s4_class(demo_peak, "GRanges")
    expect_true(length(demo_peak) > 0)
    expect_true(ncol(mcols(demo_peak)) >= 1)
})


test_that("peakAnno", {
    data("peakAnno", package = "epiSeeker")

    expect_s4_class(peakAnno, "csAnno")
    expect_s4_class(peakAnno@anno, "GRanges")
    expect_true(peakAnno@peakNum == length(peakAnno@anno))
})


test_that("peakAnnoList", {
    data("peakAnnoList", package = "epiSeeker")

    expect_type(peakAnnoList, "list")
    expect_true(length(peakAnnoList) > 0)

    lapply(peakAnnoList, function(obj) {
        expect_s4_class(obj, "csAnno")
    })
})


test_that("tagMatrix", {
    data("tagMatrix", package = "epiSeeker")

    expect_true(is.matrix(tagMatrix))
    expect_true(nrow(tagMatrix) > 0)
    expect_true(ncol(tagMatrix) > 0)
})


test_that("pwm_obj", {
    skip_if_not_installed("TFBSTools")
    data("pwm_obj", package = "epiSeeker")

    expect_true(
        inherits(pwm_obj, "PFMatrixList") ||
            inherits(pwm_obj, "PWMatrixList")
    )
    expect_true(length(pwm_obj) > 0)
})


test_that("demo_bmdata", {
    data("demo_bmdata", package = "epiSeeker")

    expect_s4_class(demo_bmdata, "bmData")

    expect_true(length(SummarizedExperiment::assays(demo_bmdata)) >= 1)
})


test_that("seq2gene_result", {
    data("seq2gene_result", package = "epiSeeker")

    expect_type(seq2gene_result, "character")
    expect_true(length(seq2gene_result) > 0)
})


test_that("gsminfo", {
    data("gsminfo", package = "epiSeeker")

    expect_true(is.data.frame(gsminfo))
    expect_true(nrow(gsminfo) > 0)

    required_cols <- c(
        "series_id", "gsm", "gpl", "organism", "title",
        "characteristics", "source_name", "extract_protocol",
        "description", "data_processing", "submission_date",
        "supplementary_file", "genomeVersion", "pubmed_id"
    )
    expect_true(all(required_cols %in% colnames(gsminfo)))
})
