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
NULL


#' Human BSseq object created by DSS
#'
#' A data set contains methylation information
#' see data-raw/Human_data_procession.R
#' @name Human_BSobj
#' @return bsseq object
NULL

#' Different methylation region created by DSS
#'
#' A data set contains information of different methylation region
#' see data-raw/Human_data_procession.R
#' @format A data frame with 29 row and 9 variables
#' \describe{
#'   \item{chr}{chromosome, the chromosome information of dmR}
#'   \item{start}{the start site of dmR}
#'   \item{end}{the end site of dmR}
#'   \item{length}{the length of dmR}
#'   \item{nCG}{Number of CpG sites contained in the DMR}
#'   \item{meanMethy1}{Average methylation levels in two conditions}
#'   \item{meanMethy2}{Average methylation levels in two conditions}
#'   \item{diff.Methy}{The difference in the methylation levels between two conditions}
#'   \item{areaStat}{	The sum of the test statistics of all CpG sites within the DMR}
#' }
#' @name Human_dmR
#' @return data frame
NULL
