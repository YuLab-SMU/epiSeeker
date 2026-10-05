library(epiSeeker)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

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

test_that("all peaks dropped reports the seqlevels mismatch", {
    ## issue #238 of ChIPseeker: no seqlevel of the peaks matches TxDb, so all
    ## peaks are dropped and the annotation code then failed with
    ## "Error: invalid subscript"
    peak <- GRanges("chr1_gl000191_random", IRanges(1000, 1200))
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

    expect_error(
        suppressWarnings(annotateSeq(peak, TxDb = txdb, verbose = FALSE)),
        "seqlevels"
    )
})

test_that("an empty input returns an empty annotation", {
    ## a zero-length peak set used to fail with "Error: invalid subscript"
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    pa <- annotateSeq(GRanges(), TxDb = txdb, verbose = FALSE)

    expect_s4_class(pa, "csAnno")
    expect_equal(length(pa@anno), 0)
    expect_equal(pa@peakNum, 0)
    expect_equal(nrow(as.data.frame(pa)), 0)
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

test_that("transcript annotation metadata follows the overlapping isoform", {
    ## issue #252: a peak can overlap one isoform while another nested
    ## transcript has the closer TSS.  All transcript-level fields must then
    ## describe the isoform that supplied the genomic annotation.
    features <- GRanges(
        "chr1",
        IRanges(start = c(100, 500), end = c(1000, 900)),
        strand = c("+", "-"),
        tx_id = c(101L, 202L),
        gene_id = c("GENE1", "GENE1")
    )
    peak <- GRanges("chr1", IRanges(700, 710), strand = "*")

    aligned <- epiSeeker:::.align_transcript_annotation(
        peak, features, index = 1L, distance = 999,
        annotationFeatureId = "202"
    )

    expect_equal(aligned$index, 2L)
    expect_equal(as.character(mcols(features)$gene_id[aligned$index]), "GENE1")
    expect_equal(as.integer(mcols(features)$tx_id[aligned$index]), 202L)
    expect_equal(aligned$distance, 190)
})

test_that("annotation and transcript-level columns describe the same transcript", {
    ## issue #252: this peak is located inside an intron of the long BRCA1
    ## isoform while the closest TSS belongs to another isoform.  The reported
    ## gene/transcript columns have to follow the transcript named in the
    ## annotation instead of the nearest one.
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    peak <- GRanges("chr1", IRanges(243832958, 243833008))

    pa <- annotateSeq(peak, TxDb = txdb, tssRegion = c(-3000, 3000),
                      level = "transcript", verbose = FALSE)
    df <- as.data.frame(pa)

    ## the annotation text names the transcript that overlaps the peak ...
    expect_match(df$annotation, "^Intron ")
    annotatedTx <- sub("^[A-Za-z']+ \\(([^/]+)/.*$", "\\1", df$annotation)

    ## ... and all other transcript-level fields come from that transcript
    expect_equal(df$transcriptId, annotatedTx)
    expect_equal(df$geneChr, as.character(seqnames(peak)))

    tx <- transcripts(txdb)
    tx <- tx[tx$tx_name == df$transcriptId]
    expect_equal(length(tx), 1L)
    expect_true(overlapsAny(peak, tx))
    expect_equal(df$geneStart, start(tx))
    expect_equal(df$geneEnd, end(tx))
    expect_equal(as.character(df$geneStrand), as.character(strand(tx)))

    ## distanceToTSS is the strand aware distance to the TSS of that transcript
    strandTx <- as.character(strand(tx))
    tss <- ifelse(strandTx == "+", start(tx), end(tx))
    dStart <- ifelse(strandTx == "+", start(peak) - tss, tss - start(peak))
    dEnd <- ifelse(strandTx == "+", end(peak) - tss, tss - end(peak))
    expect_equal(df$distanceToTSS, ifelse(abs(dStart) <= abs(dEnd),
                                          dStart, dEnd))
})

test_that("gene level annotation is aligned by gene_id", {
    ## issue #252 (gene level): .align_annotation_feature() also aligns the
    ## nearest gene with the gene of an exon/intron hit
    features <- GRanges(
        "chr1",
        IRanges(start = c(100, 500), end = c(1000, 900)),
        strand = c("+", "-"),
        gene_id = c("GENE1", "GENE2")
    )
    peak <- GRanges("chr1", IRanges(700, 710), strand = "*")

    aligned <- epiSeeker:::.align_annotation_feature(
        peak, features, index = 1L, distance = 999,
        annotationFeatureId = "GENE2", idColumn = "gene_id"
    )

    expect_equal(aligned$index, 2L)
    expect_equal(aligned$distance, 190)

    ## unknown ids leave the nearest feature untouched
    unchanged <- epiSeeker:::.align_annotation_feature(
        peak, features, index = 1L, distance = 999,
        annotationFeatureId = NA_character_, idColumn = "gene_id"
    )
    expect_equal(unchanged$index, 1L)
    expect_equal(unchanged$distance, 999)
})

test_that("gene-level annotation and geneId describe the same gene", {
    ## issue #252 (gene level): the annotation comes from an exon/intron hit
    ## while the nearest gene is found independently, so both have to follow
    ## the gene of the annotation
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    g <- suppressMessages(genes(txdb))
    gs <- g[as.character(seqnames(g)) == "chr17" & width(g) > 20000]
    strandG <- as.character(strand(gs))
    tss <- ifelse(strandG == "+", start(gs), end(gs))
    pos <- round(tss + rep(c(0.15, 0.4, 0.6, 0.85), each = length(gs)) *
                     width(gs) * ifelse(strandG == "+", 1, -1))
    own <- ifelse(strandG == "+", 1, -1) * (pos - tss)

    d <- abs(outer(pos, tss, "-"))
    diag(d) <- Inf
    cand <- unique(pos[apply(d, 1, min) < abs(own)])
    expect_true(length(cand) > 0)

    peaks <- GRanges("chr17", IRanges(utils::head(cand, 20), width = 200))
    pa <- annotateSeq(peaks, TxDb = txdb, tssRegion = c(-3000, 3000),
                      level = "gene", verbose = FALSE)
    df <- as.data.frame(pa)

    m <- match(as.character(df$geneId), names(g))
    expect_false(any(is.na(m)))
    expect_equal(as.character(df$geneId), as.character(names(g)[m]))
    expect_equal(df$geneChr, as.character(seqnames(g))[m])
    expect_equal(df$geneStart, start(g)[m])
    expect_equal(df$geneEnd, end(g)[m])
    expect_equal(as.character(df$geneStrand), as.character(strand(g))[m])

    strandRep <- as.character(strand(g))[m]
    tssRep <- as.numeric(ifelse(strandRep == "+", start(g)[m], end(g)[m]))
    dStart <- ifelse(strandRep == "+", start(peaks) - tssRep,
                     tssRep - start(peaks))
    dEnd <- ifelse(strandRep == "+", end(peaks) - tssRep, tssRep - end(peaks))
    expect_equal(df$distanceToTSS, ifelse(abs(dStart) <= abs(dEnd),
                                          dStart, dEnd))

    genic <- grepl("^(Exon|Intron)", df$annotation)
    expect_true(any(genic))
    expect_true(all(overlapsAny(peaks[genic], g[m[genic]])))
    annGene <- sub("^[A-Za-z']+ \\([^/]+/([^,]+),.*$", "\\1",
                   df$annotation[genic])
    expect_equal(annGene, as.character(df$geneId[genic]))
})

