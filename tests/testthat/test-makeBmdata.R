library(epiSeeker)
library(BSgenome.Hsapiens.UCSC.hg38)

context("test function for base modification Data")

test_that("makeBmDataFromData works for GRanges", {
    gr <- GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges(start = c(10, 20), end = c(10, 20)),
        score = c(5, 10)
    )

    bm <- makeBmDataFromData(gr, sampleNames = "s1")

    expect_s4_class(bm, "bmData")
})

test_that("makeBmDataFromData works for data.frame", {
    df <- data.frame(
        chr = c("chr1", "chr1"),
        pos = c(10, 20),
        val = c(1, 2)
    )

    bm <- makeBmDataFromData(df, sampleNames = "sampleX")

    expect_s4_class(bm, "bmData")
})


test_that("makeBmDataFromData works for list of data.frame", {
    df1 <- data.frame(chr = "chr1", pos = c(10, 20), v1 = c(1, 2))
    df2 <- data.frame(chr = "chr1", pos = c(10, 20), v1 = c(3, 4))
    lst <- list(df1, df2)

    bm <- makeBmDataFromData(lst, sampleNames = c("S1", "S2"))

    expect_s4_class(bm, "bmData")
})


test_that("plotBmProf works with list input", {
    data("demo_bmdata", package = "epiSeeker")

    bm_df <- getBmMatrix(
        region = data.frame(chr = "chr22", start = 10525991, end = 10526342),
        input = demo_bmdata,
        BSgenome = BSgenome.Hsapiens.UCSC.hg38,
        base = "C",
        motif = "CG",
        cover_depth = TRUE
    )

    df_list <- list(a = bm_df, b = bm_df)


    p <- plotBmProf(
        df = df_list,
        interactive = FALSE,
        ncol = 1
    )


    expect_true(
        inherits(p, "ggplot")
    )
})
