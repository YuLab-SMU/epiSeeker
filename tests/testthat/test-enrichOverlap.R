library(epiSeeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

context("test function for enrichAnnoOverlap ")

test_that("enrichAnnoOverlap works with example peakfile", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    peakfile <- system.file("extdata", "demo_peak.txt", package = "epiSeeker")

    res <- enrichAnnoOverlap(peakfile, peakfile, txdb)

    expect_s3_class(res, "data.frame")
    expect_true(all(c("qSample", "tSample", "qLen", "tLen", "N_OL", "pvalue", "p.adjust") %in% names(res)))
    expect_equal(nrow(res), 1)
})


test_that("enrichPeakOverlap works with GRanges and file input", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    peakfile <- system.file("extdata", "demo_peak.txt", package = "epiSeeker")
    gr <- readPeakFile(peakfile)[1:10]

    # reduce shuffle = 5 to speed up test
    res <- enrichPeakOverlap(gr, peakfile, txdb, nShuffle = 5, mc.cores = 1, verbose = FALSE)

    expect_s3_class(res, "data.frame")
    expect_true(all(c("qSample", "tSample", "qLen", "tLen", "N_OL", "pvalue", "p.adjust") %in% names(res)))
    expect_equal(nrow(res), 1)
})


test_that("shuffle", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    p <- GRanges(
        seqnames = c("chr1", "chr3"),
        ranges = IRanges(
            start = c(1, 100),
            end = c(50, 130)
        )
    )

    res <- shuffle(p, TxDb = txdb)

    expect_s4_class(res, "GRanges")
    expect_equal(length(res), 2)
    expect_true(all(width(res) == width(p)))
})


test_that("enrichOverlap.peak.internal", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    query <- GRanges("chr1", IRanges(c(100, 200), c(150, 250)))
    target <- list(
        GRanges("chr1", IRanges(120, 160)),
        GRanges("chr1", IRanges(1000, 1100))
    )

    res <- epiSeeker:::enrichOverlap.peak.internal(
        query.gr = query,
        target.gr = target,
        TxDb = txdb,
        nShuffle = 5,
        mc.cores = 1,
        verbose = FALSE
    )

    expect_true(is.list(res))
    expect_true(all(c("pvalue", "overlap") %in% names(res)))
    expect_length(res$pvalue, length(target))
})
