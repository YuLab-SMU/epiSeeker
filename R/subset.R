#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @method subset csAnno
#' @export
subset.csAnno <- function(x, subset, ...) {
    if (missing(subset)) {
        keep <- rep(TRUE, length(x@anno))
    } else {
        anno_df <- as.data.frame(x@anno)
        keep <- eval(substitute(subset), anno_df, parent.frame())

        if (!is.logical(keep)) {
            stop("'subset' must evaluate to a logical vector")
        }
        if (length(keep) == 1L) {
            keep <- rep(keep, length(x@anno))
        }
        if (length(keep) != length(x@anno)) {
            stop("'subset' must have length 1 or the same length as x@anno")
        }
        keep[is.na(keep)] <- FALSE
    }

    x@anno <- x@anno[keep]

    if (nrow(x@detailGenomicAnnotation) == length(keep)) {
        x@detailGenomicAnnotation <- x@detailGenomicAnnotation[keep, , drop = FALSE]
    }

    x@annoStat <- getGenomicAnnoStat(x@anno)
    x@peakNum <- length(x@anno)

    return(x)
}
