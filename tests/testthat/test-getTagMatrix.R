library(epiSeeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

context("test function for getTagMatrix")

test_that("getTagMatrix function for single peak file",{
  
  data(demo_peak)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  # make window through txdb object
  tagMatrix <- getTagMatrix(demo_peak, type = "start_site", by = "gene", 
                            upstream = 3000, downstream = 3000,
                            TxDb = txdb, weightCol = "V7")
  
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
