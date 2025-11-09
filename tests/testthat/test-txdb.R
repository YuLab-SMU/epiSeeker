library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(epiSeeker)
library(yulab.utils)

context("TXDB")

test_that("txdb", {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    epiSeeker:::.epiSeekerEnv(txdb)
    # test env is hg38 or not
    expect_equal(epiSeeker:::get_env_genome(), "hg38")
    expect_equal(epiSeeker:::IDType(txdb), "Entrez Gene ID")

    if(packageVersion("TxDb.Hsapiens.UCSC.hg38.knownGene") >= "3.22"){
        expect_equal(epiSeeker:::TXID2EG("70455"), "ENST00000479182.5/91653")
        expect_equal(epiSeeker:::TXID2EG("70455", geneIdOnly=TRUE), "91653")
    }
})

