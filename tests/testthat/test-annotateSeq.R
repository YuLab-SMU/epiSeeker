library(epiSeeker)
library(GenomicRanges)

context("annotateSeq metadata columns and dropped peaks")

## 用自定义 GRanges 作 TxDb，无需额外的 TxDb 数据包

test_that("geneChr / geneStrand are characters, not factor codes", {
    ## as.data.frame() turns 'seqnames' / 'strand' into factors and assigning a
    ## factor into mcols() dropped the class, so geneChr / geneStrand came out
    ## as integers (e.g. 1/2 instead of chr1/chr2 and +/-)
    peak <- GRanges(c("chr1", "chr2"),
                    IRanges(start = c(1500, 1500), width = 100),
                    strand = "+")

    txdb_mock <- GRanges(c("chr1", "chr2"),
                         IRanges(start = c(1000, 1000), width = 1000),
                         strand  = c("-", "+"),
                         gene_id = c("GENE1", "GENE2"))

    res <- annotateSeq(peak, TxDb = txdb_mock, verbose = FALSE)
    m <- S4Vectors::mcols(res@anno)

    expect_true(is.character(m$geneChr))
    expect_true(is.character(m$geneStrand))
    expect_equal(as.character(m$geneChr), c("chr1", "chr2"))
    expect_equal(as.character(m$geneStrand), c("-", "+"))
})

test_that("peaks without any feature in TxDb are reported, not silently dropped", {
    peak <- GRanges(c("chr1", "chrM"),
                    IRanges(start = c(1500, 1500), width = 100))

    txdb_mock <- GRanges("chr1", IRanges(start = 1000, end = 2000))

    expect_warning(
        annotateSeq(peak, TxDb = txdb_mock, verbose = FALSE),
        "peaks were dropped"
    )
})
