#' vennpie method generics
#'
#'
#' @docType methods
#' @name vennpie
#' @rdname vennpie-methods
#' @examples 
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' peakfile <- system.file("extdata", "sample_peaks.txt", package="epiSeeker")
#' peakAnno <- annotateSeq(peakfile, tssRegion=c(-3000, 3000), TxDb=txdb)
#' vennpie(peakAnno)
#' @export
setGeneric("vennpie", 
  function(x, r = 0.2, cex = 1.2, ...) 
  standardGeneric("vennpie")
)


#' plotDistToTSS method generics
#'
#'
#' @docType methods
#' @name plotDistToTSS
#' @rdname plotDistToTSS-methods
#' @examples
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' peakfile <- system.file("extdata", "sample_peaks.txt", package="epiSeeker")
#' peakAnno <- annotateSeq(peakfile, TxDb=txdb)
#' plotDistToTSS(peakAnno)
#' @export
setGeneric("plotDistToTSS", 
  function(x, 
    distanceColumn="distanceToTSS",
    xlab="", ylab="Binding sites (%) (5'->3')",
    title="Distribution of transcription factor-binding loci relative to TSS", 
    ...)
  standardGeneric("plotDistToTSS")
)

#' plotAnnoBar method generics
#'
#'
#' @docType methods
#' @name plotAnnoBar
#' @rdname plotAnnoBar-methods
#' @examples
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' peakfile <- system.file("extdata", "sample_peaks.txt", package="epiSeeker")
#' peakAnno <- annotateSeq(peakfile, TxDb=txdb)
#' plotAnnoBar(peakAnno)
#' @export
setGeneric("plotAnnoBar", 
  function(x,
    xlab="",
    ylab="Percentage(%)",
    title="Feature Distribution",
    ...)
  standardGeneric("plotAnnoBar")
)


#' plotAnnoPie method generics
#'
#'
#' @docType methods
#' @name plotAnnoPie
#' @rdname plotAnnoPie-methods
#' @export
setGeneric("plotAnnoPie", 
  function(x, 
    ndigit=2,
    cex=0.9,
    col=NA,
    legend.position="rightside",
    pie3D=FALSE,
    radius=0.8,
    ...)
  standardGeneric("plotAnnoPie")
)


#' getBmMatrix methods generics
#'
#'
#' @docType methods
#' @name getBmMatrix
#' @rdname getBmMatrix-methods
#' @importFrom methods setGeneric
#' @examples 
#' data(Human_BSobj)
#' require(BSgenome.Hsapiens.UCSC.hg19)
#' BSgenome_hg19 <- BSgenome.Hsapiens.UCSC.hg19
#' bmMatrix <- getBmMatrix(region = data.frame(chr = "chr1", start = 894849, end = 895849),
#'                         BSgenome = BSgenome_hg19,
#'                         input = Human_BSobj[,c(1)],
#'                         base = "C",
#'                         motif = c("CG","CHH","CHG"))
#' @export
setGeneric("getBmMatrix",
           function(region,
                    input,
                    BSgenome,
                    base = NULL,
                    motif = NULL,
                    position_bias = NULL,
                    ...){

             standardGeneric("getBmMatrix")

           })

#' makeBmDataFromData method generics
#'
#'
#' @docType methods
#' @name makeBmDataFromData
#' @rdname makeBmDataFromData-methods
#' @importFrom methods setGeneric
#' @return bmData
#' @export
setGeneric("makeBmDataFromData", function(data,
                                          sampleNames=NULL){
  standardGeneric("makeBmDataFromData")
})
