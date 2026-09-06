library(epiSeeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

context("test function for getTagMatrix")

test_that("getTagMatrix function for single peak file", {
    data(demo_peak)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    # make window through txdb object
    tagMatrix <- getTagMatrix(demo_peak,
        type = "start_site", by = "gene",
        upstream = 3000, downstream = 3000,
        TxDb = txdb, weightCol = "V7"
    )

    # test tagmatrix record
    expect_equal(attr(tagMatrix, "type"), "start_site")
    expect_equal(attr(tagMatrix, "by"), "gene")

    # test tagMatrix upstream and downstream
    expect_equal(dim(tagMatrix)[2], 6001)
})

# test the getPromoters can run normally or not
test_that("getPromoters runs with default parameters", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    res <- getPromoters(txdb)

    expect_s4_class(res, "GRanges")
    expect_true(length(res) > 0)
})

make_issue236_case <- function(type = c("site", "body")) {
    type <- match.arg(type)
    peak_seqlevels <- c("chr1", "chr2", "chr10")
    peak_starts <- c(101L, 201L, 301L)

    if (type == "site") {
        peak <- GenomicRanges::GRanges(
            seqnames = S4Vectors::Rle(factor(
                peak_seqlevels,
                levels = peak_seqlevels
            )),
            ranges = IRanges::IRanges(peak_starts, width = 5L),
            strand = "*"
        )
        S4Vectors::mcols(peak)$weight <- c(1, 2, 10)

        window_index <- c(1L, 3L, 2L)
        window_seqlevels <- peak_seqlevels[window_index]
        windows <- GenomicRanges::GRanges(
            seqnames = S4Vectors::Rle(factor(
                window_seqlevels,
                levels = window_seqlevels
            )),
            ranges = IRanges::IRanges(
                peak_starts[window_index],
                width = 5L
            ),
            strand = "+"
        )
        attr(windows, "type") <- "start_site"
        attr(windows, "upstream") <- 2
        attr(windows, "downstream") <- 2
        attr(windows, "label") <- "center"
    } else {
        peak <- GenomicRanges::GRanges(
            seqnames = S4Vectors::Rle(factor(
                peak_seqlevels,
                levels = peak_seqlevels
            )),
            ranges = IRanges::IRanges(peak_starts, width = 20L),
            strand = "*"
        )
        S4Vectors::mcols(peak)$weight <- c(1, 2, 10)

        window_index <- rep(c(1L, 3L, 2L), each = 2L)
        window_seqlevels <- peak_seqlevels[window_index]
        windows <- GenomicRanges::GRanges(
            seqnames = S4Vectors::Rle(factor(
                window_seqlevels,
                levels = unique(window_seqlevels)
            )),
            ranges = IRanges::IRanges(
                peak_starts[window_index] + rep(c(0L, 10L), 3L),
                width = 10L
            ),
            strand = "+"
        )
        attr(windows, "type") <- "body"
        attr(windows, "upstream") <- 0
        attr(windows, "downstream") <- 0
        attr(windows, "label") <- c("SS", "TS")
    }

    attr(windows, "by") <- "issue236_demo"
    list(peak = peak, windows = windows)
}

test_that("getTagMatrix keeps site windows aligned across seqlevel orders", {
    case <- make_issue236_case("site")

    tag_matrix <- getTagMatrix(
        peak = case$peak,
        windows = case$windows,
        weightCol = "weight",
        verbose = FALSE,
        ignore_strand = TRUE
    )

    expect_equal(unname(rowSums(tag_matrix)), c(5, 50, 10))
})

test_that("getTagMatrix keeps body windows aligned across seqlevel orders", {
    case <- make_issue236_case("body")

    tag_matrix <- getTagMatrix(
        peak = case$peak,
        windows = case$windows,
        weightCol = "weight",
        nbin = 5,
        verbose = FALSE,
        ignore_strand = TRUE
    )

    expect_equal(
        unname(rowSums(tag_matrix)),
        c(5, 5, 50, 50, 10, 10)
    )
})
