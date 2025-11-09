#' vennpie method generics
#'
#'
#' @docType methods
#' @name vennpie
#' @rdname vennpie-methods
#' @examples 
#' data(peakAnno)
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
#' data(peakAnno)
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
#' data(peakAnno)
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
#' require(BSgenome.Hsapiens.UCSC.hg38)
#' data(demo_bmdata)
#' bmMatrix <- getBmMatrix(region = data.frame(chr = "chr22", start = 10525991, end = 10526342),
#'                         BSgenome = BSgenome.Hsapiens.UCSC.hg38,
#'                         input = demo_bmdata,
#'                         base = "C",
#'                         motif = c("CG"))
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
