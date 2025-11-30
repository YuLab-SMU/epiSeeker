library(epiSeeker)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

context("test function for annotation")

test_that("annotateSeq handles TxDb in grange format", {

    peak <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges   = IRanges::IRanges(start = 1000, end = 1100),
        strand   = "+"
    )

    txdb_mock <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges   = IRanges::IRanges(start = 500, end = 2000),
        strand   = "+",
        gene_id  = "GENE1"
    )


    res <- annotateSeq(peak, TxDb = txdb_mock, verbose = FALSE)


    expect_s4_class(res, "csAnno")
    expect_true(length(res@anno) == 1)
    expect_equal(res@peakNum, 1)

    expect_false(res@hasGenomicAnnotation)
})



test_that("annotateSeq switches to USER_DEFINED of  TxDb in grange format", {

    peak <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges   = IRanges::IRanges(start = 1500, end = 1600)
    )

    txdb_mock <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges   = IRanges::IRanges(start = 1000, end = 2000)
    )

    res <- annotateSeq(peak, TxDb = txdb_mock, verbose = FALSE)

    expect_s4_class(res, "csAnno")
    expect_equal(res@level, "USER_DEFINED")
    expect_false(res@hasGenomicAnnotation)
    expect_equal(res@peakNum, length(peak))
    expect_false("annotation" %in% names(S4Vectors::mcols(res@anno)))
})



test_that("dropAnno by distance cutoff", {

  peak1 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = 1000, end = 1100),
    distanceToTSS = 500
  )

  peak2 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = 50000, end = 50100),
    distanceToTSS = 20000   # too far
  )

  anno_mock <- c(peak1, peak2)
  mcols(anno_mock)$distanceToTSS <- c(500, 20000)

  # Construct csAnno object
  cs <- new("csAnno",
            anno = anno_mock,
            tssRegion = c(-3000, 3000),
            level = "gene",
            hasGenomicAnnotation = FALSE,
            peakNum = 2)

  # ---- dropAnno ----
  cs2 <- dropAnno(cs, distanceToTSS_cutoff = 10000)

  expect_s4_class(cs2, "csAnno")
  expect_equal(cs2@peakNum, 1)
  expect_true(length(cs2@anno) == 1)
  expect_equal(mcols(cs2@anno)$distanceToTSS, 500)
})
