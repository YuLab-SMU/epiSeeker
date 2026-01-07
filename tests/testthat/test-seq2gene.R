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
