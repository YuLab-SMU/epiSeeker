#' @title extend filter to Peak (GRanges class object)
#' @method filter GRanges
#' @param .data granges object
#' @param ... additional parameters
#' @param .by Optional grouping variable(s) (column name or variable expression) 
#'   specifying which columns to group by when applying filters
#' @param .preserve Logical value indicating whether to preserve the original grouping 
#'   structure when .by is specified. If TRUE, group order and identities are maintained
#' @return A filtered GRanges object containing only rows that meet the specified criteria
#' @importFrom dplyr filter
#' @examples 
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' peakfile <- system.file("extdata", "sample_peaks.txt", package="epiSeeker")
#' peak <- readPeakFile(peakfile)
#' dplyr::filter(peak, fold_enrichment > 20)
#' @return grange object
#' @export
filter.GRanges <- function(.data, ..., .by = NULL, .preserve = FALSE) {
  dots = rlang::quos(...)
  as.data.frame(.data) |> 
    dplyr::filter(!!!dots, .by = .by, .preserve = .preserve) |> 
    droplevels() |> 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}

#' @title extend mutate to Peak (GRanges class object)
#' @method mutate GRanges
#' @param .data granges object
#' @param ... additional parameters
#' @param .by Optional grouping variable(s) (column name or variable expression) 
#'   specifying which columns to group by for operations
#' @param .keep Character vector specifying which columns to retain. Possible values:
#'   "all" (retain all columns, default), "used" (retain only columns used in calculations),
#'   "unused" (retain only columns not used in calculations), "none" (retain only newly created columns)
#' @param .before Column name or position index specifying where to insert new columns before
#' @param .after Column name or position index specifying where to insert new columns after
#' @return A processed GRanges object containing the added or modified columns
#' @importFrom dplyr mutate
#' @examples 
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' peakfile <- system.file("extdata", "sample_peaks.txt", package="epiSeeker")
#' peak <- readPeakFile(peakfile)
#' dplyr::mutate(peak, score = tags)
#' @export
mutate.GRanges <- function(.data, ..., .by = NULL, 
                           .keep = c("all", "used", "unused", "none"),
                           .before = NULL,
                           .after = NULL) {
  dots = rlang::quos(...)
  df = as.data.frame(.data)
  
  if (!is.null(.before) && !is.null(.after)) {
    stop("You can't supply both `.before` and `.after`.")
  }
  
  if (!is.null(.before)) {
    df = df |> 
      dplyr::mutate(!!!dots, .by = .by, .keep = .keep, .before = .before)
  } else if (!is.null(.after)) {
    df = df |> 
      dplyr::mutate(!!!dots, .by = .by, .keep = .keep, .after = .after)
  } else {
    df = df |> dplyr::mutate(!!!dots, .by = .by, .keep = .keep)
  }
  
  df |> 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}

#' @title Rename columns of a GRanges object
#' @method rename GRanges
#' @param .data A GRanges object.
#' @param ... Rename expressions in the form new_name = old_name.
#' @importFrom rlang quos
#' @importFrom dplyr rename
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @examples 
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' peakfile <- system.file("extdata", "sample_peaks.txt", package="epiSeeker")
#' peak <- readPeakFile(peakfile)
#' dplyr::rename(peak, tag = tags)
#' @return A GRanges object with renamed metadata columns.
#' @export
rename.GRanges <- function(.data, ...){
  dots <- rlang::quos(...)
  new <- as.data.frame(.data) |> 
          dplyr::rename(!!!dots) |> 
          GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  return(new)
}

#' @title arrange granges object
#' @method arrange GRanges
#' @importFrom dplyr arrange
#' @param .data granges object
#' @param ... additional parameters
#' @param .by_group If TRUE, will sort first by grouping variable. Applies to grouped data frames only.
#' @examples 
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' peakfile <- system.file("extdata", "sample_peaks.txt", package="epiSeeker")
#' peak <- readPeakFile(peakfile)
#' dplyr::arrange(peak, seqnames)
#' @return grange object
#' @export
arrange.GRanges <- function(.data, ..., .by_group = FALSE){
  dots = rlang::quos(...)
  as.data.frame(.data) |> 
    dplyr::arrange(!!!dots, .by_group = .by_group) |> 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}


