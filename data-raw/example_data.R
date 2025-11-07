# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# peakfile <- system.file("extdata", "sample_peaks.txt", package="epiSeeker")
# peakAnno <- annotateSeq(peakfile, TxDb=txdb)

# usethis::use_data(peakAnno, overwrite = TRUE, compress = "xz")


library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
set.seed(929)
# data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6418464
data <- readPeakFile("./GSM6418464_H441_FOXA2_Rep1.narrowPeak.bed.gz")
demo_peak <- data[sample(length(data), 50)]
peakAnno <- annotateSeq(demo_peak, TxDb=txdb)

write.table(demo_peak, file = "./inst/extdata/demo_peak.txt", quote = FALSE, sep = "\t")

usethis::use_data(peakAnno, overwrite = TRUE, compress = "xz")
usethis::use_data(demo_peak, overwrite = TRUE, compress = "xz")

library(RSQLite)
library(TFBSTools)
opts_base <- list()
opts_base[["collection"]] <- "CORE"
opts_base[["all_versions"]] <- FALSE
opts_base[["species"]] <- "Drosophila melanogaster"
opts_base[["tax_group"]] <- "insects"
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), JASPAR2024::db(JASPAR2024::JASPAR2024()))
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
files <- getSampleFiles()
peakAnnoList <- lapply(files, annotateSeq, TxDb=txdb, annoDb="org.Hs.eg.db",
                       tssRegion=c(-3000, 3000), verbose=FALSE)
names(peakAnnoList) <- names(peakfiles)
usethis::use_data(peakAnnoList, overwrite = TRUE, compress = "xz")
