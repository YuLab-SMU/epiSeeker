library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(epiSeeker)
library(yulab.utils)

context("TXDB")

test_that("Update txdb", {

    hg19_txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    epiSeeker:::.epiSeekerEnv(hg19_txdb)
    expect_equal(epiSeeker:::get_env_genome(), "hg19")

    hg38_txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    epiSeeker:::.epiSeekerEnv(hg38_txdb)
    expect_equal(epiSeeker:::get_env_genome(), "hg38")

})

test_that("txdb", {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    epiSeeker:::.epiSeekerEnv(txdb)
    expect_equal(epiSeeker:::IDType(txdb), "Entrez Gene ID")
    expect_equal(epiSeeker:::TXID2EG("70455"), "uc002qsd.4/1")
    expect_equal(epiSeeker:::TXID2EG("70455", geneIdOnly=TRUE), "1")
})

