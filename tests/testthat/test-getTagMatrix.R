library(epiSeeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

context("test getTagMatrix() and related functions")

test_that("getTagMatrix function for single peak file",{
  
  peak <- readPeakFile(getSampleFiles()[[4]])[1:50]
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  # make window through txdb object
  tagMatrix <- getTagMatrix(peak, type = "start_site", by = "gene", 
                            upstream = 500, downstream = 500,
                            TxDb = txdb, weightCol = "V5")
  
  expect_is(tagMatrix, "matrix")
})
