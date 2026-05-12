library(epiSeeker)

context("test function for plot venn-like fig")

test_that("vennpie.csAnno runs", {
    data(peakAnno)
    expect_s4_class(peakAnno, "csAnno")
    expect_true(peakAnno@hasGenomicAnnotation)

    expect_silent(vennpie.csAnno(peakAnno))
})


test_that("vennpie.csAnno errors on non-csAnno input", {
    fake_input <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))

    expect_error(vennpie.csAnno(fake_input))
})


test_that("vennplot works with gene list example", {
    skip_if_not_installed("ggVennDiagram")
    data(peakAnnoList, package = "epiSeeker")

    genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)

    unname_list <- unname(genes)

    expect_warning(
        p1 <- vennplot(unname_list),
        regexp = "input is not a named list"
    )

    expect_s3_class(p1, "ggplot")
})


test_that("vennplot works with named list", {
    skip_if_not_installed("ggVennDiagram")
    data(peakAnnoList, package = "epiSeeker")
    genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)

    names(genes) <- paste0("Sample", seq_along(genes))

    p2 <- vennplot(genes)

    expect_s3_class(p2, "ggplot")
})


test_that("vennplot.peakfile works with example peak files", {
    skip_if_not_installed("ggVennDiagram")
    files <- list(
        system.file("extdata", "sample_peaks.txt", package = "epiSeeker"),
        system.file("extdata", "sample_peaks.txt", package = "epiSeeker")
    )

    p3 <- vennplot.peakfile(files)

    expect_s3_class(p3, "ggplot")
})


test_that("vennplot.peakfile with custom labels", {
    skip_if_not_installed("ggVennDiagram")
    files <- list(
        system.file("extdata", "sample_peaks.txt", package = "epiSeeker"),
        system.file("extdata", "sample_peaks.txt", package = "epiSeeker")
    )

    labels <- c("A", "B")

    p4 <- vennplot.peakfile(files, labels = labels)

    expect_s3_class(p4, "ggplot")
})
