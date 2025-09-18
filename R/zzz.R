##' @importFrom yulab.utils yulab_msg
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(yulab_msg(pkgname))

  options(epiSeeker.downstreamDistance = 300)
  options(epiSeeker.ignore_1st_exon = FALSE)
  options(epiSeeker.ignore_1st_intron = FALSE)
  options(epiSeeker.ignore_downstream = FALSE)
  options(epiSeeker.ignore_promoter_subcategory= FALSE)
  
  options(aplot_align = 'y')

}

