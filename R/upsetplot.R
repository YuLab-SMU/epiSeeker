#' @importFrom grid viewport
#' @importFrom grid pushViewport
#' @importFrom grid popViewport
#' @importFrom graphics plot.new
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_minimal
#' @author Guangchuang Yu
upsetplot.csAnno <- function(x, order_by = "freq", vennpie=FALSE, vp = list(x=.6, y=.7, width=.8, height=.8)) {
    y <- x@detailGenomicAnnotation
    nn <- names(y)
    y <- as.matrix(y)

    res <- tibble::tibble(anno = lapply(seq_len(nrow(y)), function(i) nn[y[i,]]))
    g <- ggplot(res, aes_(x = ~anno)) + geom_bar() +
        xlab(NULL) + ylab(NULL) + theme_minimal() +
        ggupset::scale_x_upset(n_intersections = 20, order_by = order_by) 

    if (!vennpie) return(g)

    f <- function() vennpie(x, cex = .9)

    p <- ggplotify::as.ggplot(f) + coord_fixed() 

    ggplotify::as.ggplot(g) +
        ggimage::geom_subview(subview = p, x = vp$x, y = vp$y, width = vp$width, height = vp$height)

}
