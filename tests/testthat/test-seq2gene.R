library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(epiSeeker)


context("test function for seq2gene ")

test_that("seq2gene runs correctly on demo_peak", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    data("demo_peak", package = "epiSeeker")

    genes <- seq2gene(
        seq = demo_peak,
        tssRegion = c(-1000, 1000),
        flankDistance = 3000,
        TxDb = txdb
    )


    expect_type(genes, "character")
    expect_true(length(genes) > 0)
})

test_that("seq2gene handles regions without exon overlap", {
    ## issue #248: exons is NA when no region overlaps an exon, and
    ## `exons$gene` used to fail with
    ## "$ operator is invalid for atomic vectors"
    ## (chr1:100000-100100 overlaps introns but no exon)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peak <- GenomicRanges::GRanges("chr1",
                                   IRanges::IRanges(100000, 100100))

    genes <- seq2gene(peak, tssRegion = c(-3000, 3000),
                      flankDistance = 5000, TxDb = txdb)

    expect_type(genes, "character")
    expect_false(any(is.na(genes)))
    ## the intron host genes are reported
    expect_true(length(genes) > 0)
})

test_that("seq2gene returns an empty vector when nothing is nearby", {
    ## a peak in the chr8 centromere without any exon/intron overlap and
    ## far away from every gene
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peak <- GenomicRanges::GRanges("chr8",
                                   IRanges::IRanges(45000000, 45000100))

    empty <- seq2gene(peak, tssRegion = c(0, 0),
                      flankDistance = 1, TxDb = txdb)
    expect_equal(empty, character(0))
})
