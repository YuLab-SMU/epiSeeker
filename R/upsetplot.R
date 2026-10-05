#' @importFrom graphics plot.new
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_minimal
#' @importFrom rlang .data
#' @author Guangchuang Yu
upsetplot.csAnno <- function(x, order_by = "freq", vennpie = FALSE, vp = list(x = .6, y = .7, width = .8, height = .8)) {
    rlang::check_installed("ggupset", reason = "For upset plot.")
    y <- x@detailGenomicAnnotation
    nn <- names(y)
    y <- as.matrix(y)

    res <- tibble::tibble(anno = lapply(seq_len(nrow(y)), function(i) nn[y[i, ]]))
    g <- ggplot(res, aes(x = .data$anno)) +
        geom_bar() +
        xlab(NULL) +
        ylab(NULL) +
        theme_minimal() +
        ggupset::scale_x_upset(n_intersections = 20, order_by = order_by)

    if (!vennpie) {
        return(g)
    }

    f <- function() vennpie(x, cex = .9)

    ## The sub-view is embedded as an annotation_custom() layer, which ggplot2
    ## >= 4.0 only supports below coord_cartesian(). The former coord_fixed()
    ## therefore had to go; no replacement is added because ggplotify rasterises
    ## the grob with the aspect ratio of the device anyway and forcing a square
    ## panel only distorted it further (measured on vennpie(): anisotropy
    ## sqrt(lambda1/lambda2) 1.43 without versus 1.50 with theme(aspect.ratio=1),
    ## against 1.41 for the undistorted base graphics drawing).
    ## Otherwise the plot could not even be drawn:
    ## "`annotation_custom()` only works with `coord_cartesian()`".
    p <- ggplotify::as.ggplot(f)

    ggplotify::as.ggplot(g) +
        ggimage::geom_subview(subview = p, x = vp$x, y = vp$y, width = vp$width, height = vp$height)
}
