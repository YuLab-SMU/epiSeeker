library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(epiSeeker)
library(yulab.utils)

context("test function for TXDB")

test_that("txdb", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    epiSeeker:::.epiSeekerEnv(txdb, force = TRUE)
    # test env is hg38 or not
    expect_equal(epiSeeker:::get_env_genome(), "hg38")
    expect_equal(epiSeeker:::IDType(txdb), "Entrez Gene ID")

    ## The tx_id -> (tx_name, gene_id) mapping is annotation data that is
    ## rebuilt for every Bioconductor release. Pinning a snapshot made this
    ## test fail on the Bioc 3.24 builders: tx 70455 changed from
    ## "ENST00000479182.5/91653" to "ENST00000798887.1/102724542". Take the
    ## ids and derive the expected values from the installed TxDb instead,
    ## covering both a transcript with a gene id and one without (where
    ## TXID2EG() falls back to the transcript name).
    tx <- GenomicFeatures::transcripts(txdb, columns = c("tx_id", "gene_id"))
    ngene <- lengths(tx$gene_id)
    with_gene <- which(ngene == 1)
    without_gene <- which(ngene == 0)
    txids <- character()
    if (length(with_gene))
        txids <- c(txids, as.character(tx$tx_id[with_gene[1]]))
    if (length(without_gene))
        txids <- c(txids, as.character(tx$tx_id[without_gene[1]]))

    map <- suppressMessages(AnnotationDbi::select(
        txdb, keys = txids, columns = c("TXNAME", "GENEID"), keytype = "TXID"
    ))
    i <- match(txids, as.character(map$TXID))
    txname <- as.character(map$TXNAME[i])
    geneid <- as.character(map$GENEID[i])
    ## TXID2EG() uses the transcript name when no gene id is available
    fallback <- is.na(geneid)
    geneid[fallback] <- txname[fallback]

    expect_equal(epiSeeker:::TXID2EG(txids),
                 paste(txname, geneid, sep = "/"))
    expect_equal(epiSeeker:::TXID2EG(txids, geneIdOnly = TRUE),
                 geneid)
})
