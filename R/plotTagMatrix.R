#' @title plot peak heatmap sub functino
#' 
#' @param tagMatrix output from getTagMatrix().
#' @param xlab xlab.
#' @param ylab ylab.
#' @param title title.
#' @param palette palette to be filled in,details see \link[ggplot2]{scale_colour_brewer}.
#' @param facet_label_text_size the size of facet label text
#' @param nrow nrow to place a number of fig.
#' @param ncol ncol to place a number of fig.
#' @return ggplot object
#' @importFrom aplot plot_list
plotPeakHeatmap_sub <- function(tagMatrix,
                            xlab = "",
                            ylab = "",
                            palette = NULL,
                            title = NULL,
                            facet_label_text_size = 12,
                            nrow = NULL,
                            ncol = NULL){
    
    if(is(tagMatrix, "list")){

      multiple_samplesFlag <- attr(tagMatrix[[1]], 'multiple_samples')
      # multiple_samplesFlag 
      # True multiple samples
      # NULL multiple set

      if(is.null(multiple_samplesFlag)){
        
        # plot for multiple sets
        p <- plotPeakHeatmap_sub.internal(tagMatrix = tagMatrix,
                                          xlab = xlab,
                                          ylab = ylab,
                                          palette = palette,
                                          title = title,
                                          facet_label_text_size = facet_label_text_size)
      }else{

        # plot for multiple samples

        nc <- length(tagMatrix)
        if ( is.null(palette) || is.na(palette) ) {
          palette <- getPalette(nc)
        } else if (length(palette) != nc) {
          palette <- rep(palette[1], nc)
        } else {
          palette <- palette
        }

        if (is.null(title) || is.na(title))
          title <- names(tagMatrix)
        if (length(xlab) != nc) {
          xlab <- rep(xlab[1], nc)
        }
        if (length(ylab) != nc) {
          ylab <- rep(ylab[1], nc)
        }
        if (length(title) != nc) {
          title <- rep(title[1], nc)
        }
        tmp <- list()
    
        for (i in seq_len(nc)) {
          tmp[[i]] <- plotPeakHeatmap_sub.internal(tagMatrix = tagMatrix[[i]], 
                                                   palette = palette[i], 
                                                   xlab = xlab[i], 
                                                   ylab = ylab[i], 
                                                   title= title[i])
        }

        if(is.null(nrow) && is.null(ncol)){nrow <- 1}

        p <- plot_list(gglist = tmp,
                       ncol = ncol,
                       nrow = nrow)
        
      }

    }else{
      p <- plotPeakHeatmap_sub.internal(tagMatrix = tagMatrix,
                                        xlab = xlab,
                                        ylab = ylab,
                                        palette = palette,
                                        title = title,
                                        facet_label_text_size = facet_label_text_size)
    }

    return(p)

}

#' @title internal function of plotPeakHeatmap
#' 
#' @param tagMatrix output from getTagMatrix().
#' @param xlab xlab.
#' @param ylab ylab.
#' @param title title.
#' @param palette palette to be filled in,details see \link[ggplot2]{scale_colour_brewer}.
#' @param facet_label_text_size the size of facet label text
#' @return ggplot object
#' @importFrom yulab.utils mat2df
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 scale_fill_distiller
#' @importFrom methods missingArg
plotPeakHeatmap_sub.internal <- function(tagMatrix,
                                         xlab = "",
                                         ylab = "",
                                         palette = NULL,
                                         title = NULL,
                                         facet_label_text_size = 12){

    listFlag <- FALSE
    if(is(tagMatrix, "list")){
        listFlag <- TRUE
    }

    if(listFlag){

        by <- character(length = length(tagMatrix))
        type <- attr(tagMatrix[[1]], 'type')
        x_label <- attr(tagMatrix[[1]], 'label')
        x_label_pos <- attr(tagMatrix[[1]], 'label_position')

        tagMatrix_df <- list()
        for(idx in seq_along(tagMatrix)){

            tmp_tagMatrix <- tagMatrix[[idx]]
            by[idx] <- attr(tmp_tagMatrix, 'by')

            tmp_tagMatrix <- t(apply(tmp_tagMatrix, 1, function(x) x/max(x)))
            ii <- order(rowSums(tmp_tagMatrix))
            tmp_tagMatrix <- tmp_tagMatrix[ii,]

            colnames(tmp_tagMatrix) <- seq_len(dim(tmp_tagMatrix)[2])
            rownames(tmp_tagMatrix) <- seq_len(dim(tmp_tagMatrix)[1])

            tmp <- mat2df(tmp_tagMatrix)
            colnames(tmp) <- c("values","sample_ID","coordinate")
            tmp$sample <- names(tagMatrix)[[idx]]

            tagMatrix_df[[idx]] <- tmp
        }

        tagMatrix_df <- do.call("rbind", tagMatrix_df)

    }else{

        type <- attr(tagMatrix, 'type')
        by <- attr(tagMatrix, 'by')
        x_label <- attr(tagMatrix, 'label')
        x_label_pos <- attr(tagMatrix, 'label_position')

        tagMatrix <- t(apply(tagMatrix, 1, function(x) x/max(x)))
        ii <- order(rowSums(tagMatrix))
        tagMatrix <- tagMatrix[ii,]
  
        colnames(tagMatrix) <- seq_len(dim(tagMatrix)[2])
        rownames(tagMatrix) <- seq_len(dim(tagMatrix)[1])

        tagMatrix_df <- mat2df(tagMatrix)
        colnames(tagMatrix_df) <- c("values","sample_ID","coordinate")
    }

    if(is.null(title)){
        if(listFlag){
            title <- paste("Heatmap of peak in", gsub("_"," ", type), "of", paste(by,collapse = " and "))
        }else{
            title <- paste("Heatmap of peak in", gsub("_"," ", type), "of",by)
        }
        
    }

    sample_ID <- coordinate <- NULL

    p <- ggplot(tagMatrix_df, aes(x = coordinate,y = sample_ID)) + 
            geom_tile(aes(fill = values)) 

    if(is.null(palette)){
        p <- p + scale_fill_gradient2(low = "#91a28c", 
                                      mid = "white",
                                      high = "#8f2c37",
                                      midpoint = 0.5,     
                                      limits = c(0, 1))
    }else{
        p <- p + scale_fill_distiller(palette = palette) 
    }

    if(listFlag){
        p <-  p + facet_grid(sample ~ .,scales = "free_y",space = "free") +
                    theme(strip.text.y = element_text(color = "black",face = "bold",
                                                      size = facet_label_text_size),
                          strip.background = element_blank())
    }


    p <- p + scale_x_continuous(breaks = x_label_pos,
                                labels = x_label,expand = c(0,0))+ 
            labs(x = xlab, y = ylab, title = title) +
            theme(axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  axis.line.y = element_blank(),
                  panel.grid=element_blank(),
                  panel.background = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  panel.border = element_rect(colour = "black",         
                                              fill = NA, 
                                              linewidth = 1.2)) +
            scale_y_continuous(expand = c(0,0))
    

    return(p)

}


#' @title plot peak profile
#' 
#' @param tagMatrix output from getTagMatrix().
#' @param xlab xlab.
#' @param ylab ylab.
#' @param title title.
#' @param conf confidence interval.
#' @param facet one of 'none', 'row' and 'column'.
#' @param free_y if TRUE, y will be scaled.
#' @param statistic_method method to do statistic. one of "mean", "median", "min", "max", "sum", "std"
#' @param missingDataAsZero set missing data as zero or not. 
#' @param ... additional parameters
#' @author G Yu; Y Yan
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 scale_color_manual
#' @return ggplot object
#' @examples 
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' peakfile <- system.file("extdata", "sample_peaks.txt", package="epiSeeker")
#' peak <- readPeakFile(peakfile)
#' mt <- getTagMatrix(peak, type = "start_site", TxDb = txdb, nbin = 100)
#' plotPeakProf(mt)
#' @export 
plotPeakProf <- function(tagMatrix,
                         xlab="Genomic Region (5'->3')",
                         ylab = "Peak Count Frequency",
                         conf = NULL,
                         title = "",
                         facet="none", 
                         free_y = TRUE,
                         statistic_method = "mean",
                         missingDataAsZero = TRUE,
                         ...){
    
    multiple_samplesFlag <- attr(tagMatrix[[1]], 'multiple_samples')

    listFlag <- FALSE
    if (is(tagMatrix, "list")) {
      if (is.null(names(tagMatrix)) ) {
        nn <- paste0("peak", seq_along(tagMatrix))
        warning("input is not a named list, set the name automatically to ", paste(nn, collapse=' '))
        names(tagMatrix) <- nn
        ## stop("tagMatrix should be a named list...")
      }
      facet <- match.arg(facet, c("none", "row", "column"))
      x_label <- attr(tagMatrix[[1]], 'label')
      x_label_pos <- attr(tagMatrix[[1]], 'label_position')
      dash_pos <- attr(tagMatrix[[1]], 'dash_pos')

      listFlag <- TRUE
    }else{
      x_label <- attr(tagMatrix, 'label')
      x_label_pos <- attr(tagMatrix, 'label_position')
      dash_pos <- attr(tagMatrix, 'dash_pos') 
    }

    ## S4Vectors change the behavior of ifelse
    ## see https://support.bioconductor.org/p/70871/
    ##
    ## conf <- ifelse(missingArg(conf), NA, conf)
    ##
    conf <- if(is.null(conf)) NA else conf

      
    if (listFlag) {

      xlim <- c(1,ncol(tagMatrix[[1]]))
      tagCount <- lapply(tagMatrix, function(x) getTagCount(x, xlim = xlim, conf = conf, 
                                                            statistic_method = statistic_method, 
                                                            missingDataAsZero = missingDataAsZero, ...))
      tagCount <- list_to_dataframe(tagCount)
      tagCount$.id <- factor(tagCount$.id, levels=names(tagMatrix))

      Upper <- Lower <- .id <- NULL
      p <- ggplot(tagCount, aes(pos, group=.id, color=.id))
      if (!(is.na(conf))) {
        p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = .id),
                            linetype = 0, alpha = 0.2)
      }
    } else {
      xlim <- c(1,ncol(tagMatrix))
      tagCount <- getTagCount(tagMatrix, xlim = xlim, conf = conf, 
                              statistic_method = statistic_method, 
                              missingDataAsZero = missingDataAsZero, ...)
      p <- ggplot(tagCount, aes(pos))
      if (!(is.na(conf))) {
        p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper),
                            linetype = 0, alpha = 0.2)
      }
    }

    value <- NULL
    p <- p + geom_line(aes(y = value)) +
              geom_vline(xintercept = dash_pos, linetype="longdash") + 
              scale_x_continuous(breaks = x_label_pos,
                                 labels = x_label,expand = c(0,0)) + 
              labs(x = xlab, y = ylab, title = title) +
              theme_bw() + 
              theme(legend.title=element_blank(),
                    plot.title = element_text(hjust = 0.5))

    if(listFlag){
     
        cols <- getCols(length(tagMatrix))
        p <- p + scale_color_manual(values=cols)
        if (facet == "row") {
          if (free_y) {
            p <- p + facet_grid(.id ~ ., scales = "free_y")
          } else {
            p <- p + facet_grid(.id ~ .)
          }
        } else if (facet == "column") {
          if (free_y) {
            p <-  p + facet_grid(. ~ .id, scales = "free_y")
          } else {
            p <-  p + facet_grid(. ~ .id)
          }
        }

        if(facet != "none") {
          p <- p + theme(legend.position="none")
        }


    }
    
    return(p)
}


#' @title plotPeakHeatmap function
#' 
#' @param tagMatrix output from getTagMatrix().
#' @param plot_prof combine prof or not. Default: TRUE
#' @param xlab xlab.
#' @param ylab ylab.
#' @param title title.
#' @param conf confidence interval.
#' @param statistic_method method to do statistic. one of "mean", "median", "min", "max", "sum", "std"
#' @param missingDataAsZero set missing data as zero or not. 
#' @param facet one of 'none', 'row' and 'column'.
#' @param free_y if TRUE, y will be scaled.
#' @param palette palette to be filled in,details see \link[ggplot2]{scale_colour_brewer}.
#' @param facet_label_text_size the size of facet label text
#' @param nrow nrow to place a number of fig.
#' @param ncol ncol to place a number of fig.
#' @param height_proportion the proportion of profiling picture and heatmap
#' @param ... additional parameters
#' @importFrom aplot insert_bottom
#' @examples 
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' peakfile <- system.file("extdata", "sample_peaks.txt", package="epiSeeker")
#' peak <- readPeakFile(peakfile)
#' mt <- getTagMatrix(peak, type = "start_site", TxDb = txdb, nbin = 100)
#' plotPeakHeatmap(mt)
#' @return ggplot object
#' @export 
plotPeakHeatmap <- function(tagMatrix,
                            plot_prof = TRUE,
                            xlab = "",
                            ylab = "",
                            palette = NULL,
                            title = NULL,
                            facet_label_text_size = 12,
                            nrow = NULL,
                            ncol = NULL,
                            conf = NULL,
                            statistic_method = "mean",
                            missingDataAsZero = TRUE,
                            facet="none", 
                            free_y = TRUE,
                            height_proportion = 4,
                            ...){
    
    if(is(tagMatrix, "list")){
      
      multiple_samplesFlag <- attr(tagMatrix[[1]], 'multiple_samples')
      # multiple_samplesFlag 
      # True multiple samples
      # NULL multiple set

      if(is.null(multiple_samplesFlag)){

        if(is.null(title))
          title <- ""

        heapmap_p <- plotPeakHeatmap_sub(tagMatrix = tagMatrix,
                                         xlab = "gene distance (bp)",
                                         ylab = ylab,
                                         palette = palette,
                                         title = "",
                                         facet_label_text_size = facet_label_text_size,
                                         nrow = nrow,
                                         ncol = ncol)

        if(plot_prof){
          prof_p <- plotPeakProf(tagMatrix = tagMatrix,
                                 xlab = xlab,
                                 ylab = "Peak Count Frequency",
                                 conf = conf,
                                 statistic_method = statistic_method, 
                                 missingDataAsZero = missingDataAsZero,
                                 title = title,
                                 facet = 'none', 
                                 free_y = free_y,
                                 ...) +
                        theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank())

          p <- prof_p %>% insert_bottom(heapmap_p,height = height_proportion)
        }else{
          p <- heapmap_p
        }
        

      }else{

        tmp <- list()
        if(is.null(title)){title <- names(tagMatrix)}

        nc <- length(tagMatrix)
        if ( is.null(palette) || is.na(palette) ) {
          palette <- getPalette(nc)
        } else if (length(palette) != nc) {
          palette <- rep(palette[1], nc)
        } else {
          palette <- palette
        }

        for(i in seq_along(tagMatrix)){
          heapmap_p <- plotPeakHeatmap_sub(tagMatrix = tagMatrix[[i]],
                                           xlab = "gene distance (bp)",
                                           ylab = ylab,
                                           palette = palette[i],
                                           title = "",
                                           facet_label_text_size = facet_label_text_size,
                                           nrow = nrow,
                                           ncol = ncol)
          
          if(plot_prof){
            prof_p <- plotPeakProf(tagMatrix = tagMatrix[[i]],
                                   xlab = xlab,
                                   ylab = "Peak Count Frequency",
                                   title = title[i],
                                   conf = conf,
                                   statistic_method = statistic_method,
                                   missingDataAsZero = missingDataAsZero,
                                   facet = facet, 
                                   free_y = free_y,
                                   ...)+
                        theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank())

            tmp[[i]] <- prof_p %>% insert_bottom(heapmap_p, height = height_proportion)
          }else{
            tmp[[i]] <- heapmap_p
          }
          
        }

        if (is.null(ncol) && is.null(nrow))
          nrow <- 1
    
        p <- plot_list(gglist = tmp,
                       ncol = ncol,
                       nrow = nrow)
      }


    }else{

      if(is.null(title)){
        title <- ""
      }
      heapmap_p <- plotPeakHeatmap_sub(tagMatrix = tagMatrix,
                                       xlab = "gene distance (bp)",
                                       ylab = ylab,
                                       palette = palette,
                                       title = "",
                                       facet_label_text_size = facet_label_text_size,
                                       nrow = nrow,
                                       ncol = ncol)

      if(plot_prof){
        prof_p <- plotPeakProf(tagMatrix = tagMatrix,
                               xlab = xlab,
                               ylab = "Peak Count Frequency",
                               conf = conf,
                               statistic_method = statistic_method,
                               missingDataAsZero = missingDataAsZero,
                               title = title,
                               facet = facet, 
                               free_y = free_y,
                               ...) +
                        theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank())
    
        p <- prof_p %>% insert_bottom(heapmap_p,height = height_proportion)
      }else{
        p <- heapmap_p
      }
      
    }

    return(p)

}