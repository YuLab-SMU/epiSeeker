#' plot the overlap of a list of object
#'
#' venn plot produced through this way has colors which can be defined by users using
#' ggplot2 grammar e.g.(scale_fill_distiller()). And users can specify any details, like digital number,
#' text size and showing percentage or not, by inputting `...` extra parameters.
#' 
#' @title vennplot
#' @param Sets a list of object, can be vector or GRanges object.
#' @param ... extra parameters using ggVennDiagram. Details see \link[ggVennDiagram]{ggVennDiagram}
#' @return venn plot that summarize the overlap of peaks
#' from different experiments or gene annotation from
#' different peak files.
#' @examples
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' peakfiles <- getSampleFiles()
#' peakAnnoList <- lapply(peakfiles, annotateSeq, TxDb = txdb)
#' names(peakAnnoList) <- names(peakfiles)
#' genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
#' vennplot(genes)
#' @export
#' @author G Yu
vennplot <- function(Sets,...) {
    if (is.null(names(Sets))) {
        nn <- paste0("Set", seq_along(Sets))
        warning("input is not a named list, set the name automatically to ", paste(nn, collapse = " "))
        names(Sets) <- nn
        ## stop("input object should be a named list...")
    }

    rlang::check_installed('ggVennDiagram', reason = 'For venn plot.')
    if (requireNamespace("ggVennDiagram")){
        p <- ggVennDiagram::ggVennDiagram(Sets, ...)
    }

    return(p)

}

#' vennplot for peak files
#'
#'
#' @title vennplot.peakfile
#' @param files peak files
#' @param labels labels for peak files
#' @return figure
#' @export
#' @examples 
#' files <- list(system.file("extdata", "sample_peaks.txt", package="epiSeeker"),
#'               system.file("extdata", "sample_peaks.txt", package="epiSeeker"))
#' vennplot.peakfile(files)
#' @author G Yu
vennplot.peakfile <- function(files, labels=NULL) {
    peak.Sets <- lapply(files, readPeakFile)
    if (is.null(labels)) {
        ## remove .xls or .bed of the file names as labels
        labels <- sub("\\.\\w+$", "", files)
    }
    names(peak.Sets) <- labels
    vennplot(peak.Sets)
}


