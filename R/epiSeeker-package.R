#' @keywords internal
"_PACKAGE"



#' Information Datasets
#' 
#' ucsc genome version, precalcuated data and gsm information
#' 
#' @name gsminfo
#' @aliases ucsc_release
#' @format Data frame record ucsc genome version, precalcuated data and gsm information
#' @docType data
#' @keywords datasets
#' @return data frame
NULL

#' Example data of peak annotation
#'
#' Peak annotation result. See data-raw/example_data.R
#' @name peakAnno
#' @return csAnno object
NULL

#' demo peak file
#'
#' Peak in Grange object. See data-raw/example_data.R
#' @name demo_peak
#' @return Grange object
NULL


#' Example data of a list of peak annotation
#'
#' A list of peak annotation result. See data-raw/example_data.R
#' @name peakAnnoList
#' @return list of csAnno object
NULL

#' Example data of tagMatrix
#'
#' tagMatrix result. See data-raw/example_data.R
#' @name tagMatrix
#' @return matrix
NULL

#' motif reference for Drosophila melanogaster
#'
#' motif reference result. See data-raw/example_data.R
#' @name pwm_obj
#' @return pwm_obj
NULL

#' demo base modification data
#'
#' demo base modification data. See data-raw/example_data.R
#' @name demo_bmdata
#' @return bmData object
NULL

#' Result of seq2gene
#'
#' See data-raw/example_data.R
#' @name seq2gene_result
#' @return vector of gene names
NULL

#' Name of the epiSeeker cache environment (internal static variable)
#' @format character vector 
epiSeekerCache <- "epiSeekerEnv"