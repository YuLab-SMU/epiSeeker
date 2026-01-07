library(epiSeeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(gggenes)

context("test function for plotGeneTrack ")

test_that("plotGeneTrack", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    p <- plotGeneTrack(
        txdb = txdb,
        chr = "chr8",
        start_pos = 126712193,
        end_pos = 126713193
    )

    expect_s3_class(p, "ggplot")
})


test_that("plotGeneTrack works with palette option", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    p <- plotGeneTrack(
        txdb = txdb,
        chr = "chr8",
        start_pos = 126712193,
        end_pos = 126713193,
        palette = "Set3"
    )

    expect_s3_class(p, "ggplot")
})


test_that("plotGeneTrack works with highlight", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    p1 <- plotGeneTrack(
        txdb = txdb,
        chr = "chr8",
        start_pos = 126712193,
        end_pos = 126713193,
        highlight = c(126712500, 126712800)
    )

    expect_s3_class(p1, "ggplot")

    p2 <- plotGeneTrack(
        txdb = txdb,
        chr = "chr8",
        start_pos = 126712193,
        end_pos = 126713193,
        highlight = list(
            c(126712300, 126712500),
            c(126712900, 126713100)
        )
    )
    expect_s3_class(p2, "ggplot")
})


test_that("plotGeneTrack works with auto_x_axis", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    p <- plotGeneTrack(
        txdb = txdb,
        chr = "chr8",
        start_pos = 126712193,
        end_pos = 126713193,
        auto_x_axis = FALSE
    )

    expect_s3_class(p, "ggplot")
})


test_that("plotGeneTrack work with select_gene by gene_id", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    genes <- GenomicFeatures::genes(txdb)
    win <- GRanges("chr8", IRanges::IRanges(126712193, 126713193))
    ids <- names(subsetByOverlaps(genes, win))

    p <- plotGeneTrack(
        txdb = txdb,
        chr = "chr8",
        start_pos = 126712193,
        end_pos = 126713193,
        select_gene = ids[1]
    )
    expect_s3_class(p, "ggplot")
})
