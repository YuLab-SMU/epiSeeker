library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakfile <- system.file("extdata", "sample_peaks.txt", package="epiSeeker")
peakAnno <- annotateSeq(peakfile, TxDb=txdb)

usethis::use_data(peakAnno, overwrite = TRUE, compress = "xz")


library(RSQLite)
library(TFBSTools)
opts_base <- list()
opts_base[["collection"]] <- "CORE"
opts_base[["all_versions"]] <- FALSE
opts_base[["species"]] <- "Drosophila melanogaster"
opts_base[["tax_group"]] <- "insects"
sq24 <- DBI::dbConnect(RSQLite::SQLite(), JASPAR2024::db(JASPAR2024::JASPAR2024()))
pwm_obj <- TFBSTools::getMatrixSet(sq24, opts_base)

usethis::use_data(pwm_obj, overwrite = TRUE, compress = "xz")


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peak <- readPeakFile(getSampleFiles()[[4]])
tagMatrix <- getTagMatrix(peak, type = "start_site", by = "gene", 
                          upstream = 3000, downstream = 3000,
                          TxDb = txdb, weightCol = "V5", nbin = 500)
usethis::use_data(tagMatrix, overwrite = TRUE, compress = "xz")


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakfiles <- getSampleFiles()
peakAnnoList <- lapply(peakfiles, annotateSeq, TxDb = txdb)
names(peakAnnoList) <- names(peakfiles)
usethis::use_data(peakAnnoList, overwrite = TRUE, compress = "xz")
