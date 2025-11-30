library(epiSeeker)

context("test function for subset csAnno")

test_that("subset.csAnno", {

    data(peakAnno)

    original <- peakAnno
    n_original <- original@peakNum


    keep_chr <- as.character(GenomicRanges::seqnames(original@anno)[1])

    res <- subset(original, seqnames == keep_chr)

    expect_s4_class(res, "csAnno")
    expect_true(res@peakNum <= n_original)


    expect_true(all(GenomicRanges::seqnames(res@anno) == keep_chr))
    expect_equal(nrow(res@detailGenomicAnnotation), length(res@anno))

    expect_true(nrow(res@annoStat) > 0)
    expect_true(sum(res@annoStat$Frequency) > 0)


    original_indices <- paste(
        GenomicRanges::seqnames(original@anno),
        BiocGenerics::start(original@anno),
        BiocGenerics::end(original@anno),
        sep = "_"
    )
    res_indices <- paste(
        GenomicRanges::seqnames(res@anno),
        BiocGenerics::start(res@anno),
        BiocGenerics::end(res@anno),
        sep = "_"
    )
    expect_true(all(res_indices %in% original_indices))

})
