library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(epiSeeker)

context("TXDB")

test_that("txdb", {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    epiSeeker:::.epiSeekerEnv(txdb)
    expect_equal(epiSeeker:::IDType(txdb), "Entrez Gene ID")
    expect_equal(epiSeeker:::TXID2EG("70455"), "uc002qsd.4/1")
    expect_equal(epiSeeker:::TXID2EG("70455", geneIdOnly=TRUE), "1")
})

