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


test_that("enrichPeakOverlap accepts a single GRanges as target", {
    ## issue #84 of ChIPseeker: a bare GRanges was not wrapped into a list, so
    ## the permutation test failed with
    ## "GRanges objects don't support [[, as.list(), lapply()"
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    set.seed(1)
    q <- GRanges("chr1", IRanges(sample(1000000:1400000, 5), width = 300))
    t <- GRanges("chr1", IRanges(sample(1000000:1400000, 5), width = 300))

    res <- enrichPeakOverlap(q, t, TxDb = txdb, nShuffle = 20,
                             mc.cores = 1, verbose = FALSE)

    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 1)
    expect_equal(res$qLen, 5)
    expect_equal(res$tLen, 5)
    expect_true(res$N_OL >= 0 && res$N_OL <= 5)
})


test_that("the overlap count is symmetric while the tested ratio is not", {
    ## issue #84: N_OL is direction free, the p-value normalises by the target
    ## size and shuffles the target, so swapping the arguments can change it
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    set.seed(11)
    hot <- sample(1000000:1400000, 3)
    A <- GRanges("chr1", IRanges(c(sample(hot, 2),
                                   sample(2000000:2400000, 18)), width = 300))
    B <- GRanges("chr1", IRanges(c(sample(hot, 3),
                                   sample(3000000:3400000, 97)), width = 300))

    set.seed(1)
    ab <- enrichPeakOverlap(A, B, TxDb = txdb, nShuffle = 100, mc.cores = 1,
                            verbose = FALSE)
    set.seed(1)
    ba <- enrichPeakOverlap(B, A, TxDb = txdb, nShuffle = 100, mc.cores = 1,
                            verbose = FALSE)

    ## the number of overlapping peaks is direction free ...
    expect_equal(ab$N_OL, ba$N_OL)
    expect_equal(ab$N_OL, length(intersect(A, B)))
    expect_equal(ab$qLen, ba$tLen)
    expect_equal(ab$tLen, ba$qLen)

    ## ... while the ratio the p-value is based on is target normalised, so it
    ## changes when the two peak sets differ in size
    expect_false(isTRUE(all.equal(ab$N_OL / ab$tLen, ba$N_OL / ba$tLen)))

    expect_true(all(c(ab$pvalue, ba$pvalue) > 0))
    expect_true(all(c(ab$pvalue, ba$pvalue) <= 1))
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
