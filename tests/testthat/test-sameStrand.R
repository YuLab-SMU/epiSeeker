library(GenomicRanges)

context("sameStrand is honoured when detecting overlaps")

## NEG_gene: 反链，TSS = 300；POS_gene: 正链，TSS = 1000
features <- GRanges(
    "chr1",
    IRanges(start = c(150, 1000), end = c(300, 1200)),
    strand  = c("-", "+"),
    gene_id = c("NEG_gene", "POS_gene"),
    tx_id   = c(1L, 2L)
)

nearestGene <- function(peak_strand, overlap = "TSS", sameStrand = FALSE,
                        start = 100, end = 200) {
    peaks <- GRanges("chr1", IRanges(start = start, end = end),
                     strand = peak_strand)
    res <- epiSeeker:::getNearestFeatureIndicesAndDistances(
        peaks, features, sameStrand = sameStrand, overlap = overlap
    )
    as.character(features$gene_id[res$index])
}

test_that("sameStrand=TRUE only matches same-strand features, overlap='all'", {
    ## peak 整个落在反链 feature 的区间内
    expect_equal(nearestGene("+", overlap = "all", sameStrand = TRUE),
                 "POS_gene")
})

test_that("sameStrand=TRUE only matches same-strand features, overlap='TSS'", {
    ## peak 正好盖住反链 feature 的 TSS 点
    expect_equal(nearestGene("+", overlap = "TSS", sameStrand = TRUE,
                             start = 290, end = 310),
                 "POS_gene")
})

test_that("peaks with ambiguous strand ('*') still match any strand", {
    expect_equal(nearestGene("*", overlap = "all", sameStrand = TRUE),
                 "NEG_gene")
})

test_that("sameStrand=FALSE keeps the previous behaviour", {
    expect_equal(nearestGene("+", overlap = "all", sameStrand = FALSE),
                 "NEG_gene")
    expect_equal(nearestGene("+", overlap = "TSS", sameStrand = FALSE),
                 "NEG_gene")
})
