#' getBmMatrix method for [bsseq::BSseq]
#'
#' @docType methods
#' @rdname getBmMatrix-methods
#' @title getBmMatrix method
#' @param region base modification region in the form of data.frame, having columns of "chr","start" and "end"
#' @param input the input data stored in [bsseq::BSseq] objects or BSseqExtra objects
#' @param BSgenome genome reference
#' @param base one of A/T/G/C/U
#' @param motif the motif(e.g C:CG/CH, A:GAGG/AGG) of the base modification
#' @param position_bias 1-base bias. e.g position_bias = 1("C" in "CHH"), position_bias = 2("A" in "GAGG")
#' @param cover_depth take the depth of cover into account or not
#' @param ... other parameters
#' @aliases getBmMatrix, BSseq-methods
#' @return data.frame
#' @importFrom methods setMethod
#' @exportMethod getBmMatrix
setMethod(
    "getBmMatrix", signature(input = "BSseq"),
    function(region,
             input,
             BSgenome,
             base = NULL,
             motif = NULL,
             position_bias = NULL,
             cover_depth = TRUE,
             ...) {
        getBmMatrix.BSseq(
            region = region,
            input = input,
            BSgenome = BSgenome,
            base = base,
            motif = motif,
            position_bias = position_bias,
            cover_depth = cover_depth
        )
    }
)
