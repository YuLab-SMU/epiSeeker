#' @title plotCov
#' @details plot peak coverage
#' 
#' @param peak peak file or GRanges object.
#' @param weightCol weight column of peak.
#' @param facet_level facet_level.
#' @param highlight a region or a list of region to highlight.
#' @param highlight_color colors of highlight rect. Default "#c6c3c3"
#' @param highlight_alpha alpha of highlight rect.
#' @param xlab xlab.
#' @param ylab ylab.
#' @param interactive produce interactive fig or not.
#' @param width_svg  width_svg
#' @param height_svg height_svg
#' @param title title.
#' @param x_text_size the size of x text.
#' @param y_text_size the size of y text.
#' @param chrs selected chromosomes to plot, all chromosomes by default.
#' @param xlim ranges to plot, default is whole chromosome.
#' @param lower lower cutoff of coverage signal.
#' @param fill_color specify the color/palette for the plot. Order matters.
#' @param facet_label_text_size the size of facet label text.
#' @param facet_var how to facet. one of c("chr~.", ".~ chr", ".~.id", ".id~.", ".id~chr", "chr~.id")
#' @param facet_scales how to scale facet data. Default: "free".
#' @param legend_position legend_position
#' @param add_cluster_tree add cluster tree for samples or not.
#' @param cluster_dist_method method for calculate cluster tree. Details see stat::dist
#' @param cluster_hclust_methond method for hclust. Details see stat::hclust
#' @param add_coaccess add co-accessibility or not 
#' @param curvature curvature.
#' @param coaccess_top_n top n co-accessibility to show, default: 3.
#' @param coaccess_cor_threshold co-access peak cor threshold.
#' @param coaccess_legend_pos the legend position of co-accessibiliy plot legend.
#' @param coaccess_legend_text_size the legend position of co-accessibiliy plot legend text size.
#' @param coaccess_legend_title_size the legend position of co-accessibiliy plot legend title size.
#' @param design the design layout of figure.
#' @return ggplot2 object
#' @import GenomeInfoDb
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_blank
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 geom_curve
#' @importFrom ggiraph geom_rect_interactive
#' @importFrom ggiraph girafe
#' @importFrom ggiraph opts_hover
#' @importFrom ggiraph opts_hover_inv
#' @importFrom ggiraph opts_tooltip
#' @importFrom ggiraph opts_sizing
#' @importFrom ggiraph opts_zoom
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats as.formula
#' @export
#' @author G Yu
plotCov <- function(peak, weightCol = NULL,
                    facet_level = NULL,
                    highlight = NULL,
                    highlight_color = "#c6c3c3",
                    highlight_alpha = 0.2,
                    xlab  = "Chromosome Size (bp)",
                    ylab  = "",
                    interactive = FALSE,
                    width_svg = 10, 
                    height_svg = 6,
                    title = "ChIP Peaks over Chromosomes",
                    x_text_size = 10,
                    y_text_size = 10,
                    facet_label_text_size = 10,
                    chrs  = NULL,
                    xlim  = NULL,
                    facet_var = NULL,
                    facet_scales = "free",
                    lower = 1,
                    fill_color = "black",
                    add_cluster_tree = FALSE,
                    cluster_dist_method = "euclidean",
                    cluster_hclust_methond = "complete",
                    legend_position = NULL,
                    add_coaccess = FALSE,
                    curvature = 0.3,
                    coaccess_top_n = NULL,
                    coaccess_cor_threshold = NULL,
                    design = NULL,
                    coaccess_legend_pos = c(0.9,0.5),
                    coaccess_legend_text_size = 10,
                    coaccess_legend_title_size = 12) {
    
    facet_var <- match.arg(facet_var, c("chr~.", ".~chr", ".~.id", ".id~.", ".id~chr", "chr~.id"))

    isList <- is.list(peak)
    if(!isList) {  # Note: don't support data.frame
        tm <- getChrCov(peak = peak, weightCol = weightCol, chrs = chrs, xlim = xlim, lower = lower)
    } else {
        ltm <- lapply(peak, getChrCov, weightCol = weightCol, chrs = chrs, xlim = xlim, lower = lower)
        if (is.null(names(ltm))) {
            nn <- paste0("peak", seq_along(ltm))
            warning("input is not a named list, set the name automatically to ", paste(nn, collapse = ' '))
            names(ltm) <- nn
        }
        tm <- dplyr::bind_rows(ltm, .id = ".id")
        chr.sorted <- sortChrName(as.character(unique(tm$chr)))
        tm$chr <- factor(tm$chr, levels = chr.sorted)

        if(!is.null(facet_level)){
            tm$.id <- factor(tm$.id, levels = facet_level, labels = facet_level)
        }
    }

    if(add_cluster_tree){
        
        if(!isList) stop("add_cluster_tree function apply to a list of peak")
        if(facet_var != ".id~.") stop("add_cluster_tree function apply to a list of peak in facet_var(\".id~.\") ")

        mat <- grange2mt(peak, weightCol)
        d <- dist(t(mat), method = cluster_dist_method)
        hc <- hclust(d, method = cluster_hclust_methond)

        rlang::check_installed(c("ape", "ggtree"))

        if(requireNamespace("ape", quietly = TRUE)) tree <- ape::as.phylo(hc)
        
        if (requireNamespace("ggtree", quietly = TRUE)){
            tree_p <- ggtree::ggtree(tree) + 
                        ggplot2::scale_x_reverse()
        }

        tree_p_data <-  as.data.frame(tree_p$data)
        tree_p_data <- tree_p_data[tree_p_data$isTip,]
        facet_level <- tree_p_data[order(tree_p_data$y,decreasing = TRUE),"label"]

        tm$.id <- factor(tm$.id, levels = facet_level, labels = facet_level)
    }
    
    chr <- start <- end <- value <- .id <- NULL
    
    if(length(tm$chr) == 0){
        p <- ggplot(data.frame(x = 1)) + geom_blank()
    } else {
        p <- ggplot(tm, aes(start, value))
        
        ## p <- p + geom_segment(aes(x=start, y=0, xend=end, yend= value))
        if (isList) {
            if (length(fill_color) == length(peak) && all(is_valid_color(fill_color))){
                cols = fill_color
            } else {
                cols = generate_colors(fill_color, n = length(peak))
            }

            if(interactive){
                p <- p + geom_rect_interactive(aes(xmin = start, ymin = 0, xmax = end, 
                                                   ymax = value, fill = .id, color = .id,
                                                   tooltip = paste0("peak start: ", start, "\n",
                                                                "peak end: ", end, "\n",
                                                                "peak value: ", round(value, 2)),
                                                   data_id = paste0(chr, "_", start, "_", end)),
                                               hover_nearest = TRUE)
            }else{
                p <- p + geom_rect(aes(xmin = start, ymin = 0, xmax = end, ymax = value, fill = .id, color = .id)) 
            }
            
            
            p <- p +
                scale_color_manual(values = cols) +
                scale_fill_manual(values = cols)

        } else {
            if(interactive){
                p <- p + geom_rect_interactive(aes(xmin = start, ymin = 0, xmax = end, ymax = value,
                                       tooltip = paste0("peak start: ", start, "\n",
                                                        "peak end: ", end, "\n",
                                                        "peak value: ", round(value, 2)), 
                                        data_id = paste0(chr, "_", start, "_", end)),
                                        fill = fill_color, color = fill_color, hover_nearest = TRUE)
            }else{
                p <- p + geom_rect(aes(xmin = start, ymin = 0, xmax = end, ymax = value), fill = fill_color, color = fill_color)
            }
            
        }
        
        # if(length(unique(tm$chr)) > 1) {
        #     p <- p + facet_grid(chr ~., scales="free")
        # }
        
    }

    if(length(unique(tm$chr)) > 1 || length(peak) > 1){
        if(is.null(facet_var)){
            if(length(unique(tm$chr)) > 1 && length(peak) > 1){
                facet_var <- "chr ~ .id"
            }else if(length(unique(tm$chr)) > 1){
                facet_var <- "chr ~."
            }else{
                facet_var <- ".id ~."
            }
        }

        p <- p + facet_grid(as.formula(facet_var), scales = facet_scales)
        

    }
    
    p <- p + theme_classic()
    p <- p + labs(x = xlab, y = ylab, title = title, fill = NULL, color = NULL)
    p <- p + scale_y_continuous(expand = c(0,0))
    p <- p + theme(strip.text.y=element_text(angle=360,
                                             size = facet_label_text_size),
                   axis.text.x = element_text(size = x_text_size),
                   axis.text.y = element_text(size = y_text_size))
    p <- p + scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_si("")))
    
    

    if(!is.null(legend_position)){
        p <- p + theme(legend.position = legend_position)
    }

    if (!is.null(xlim) && !all(is.na(xlim)) && is.numeric(xlim) && length(xlim) == 2) {
        p <- p + xlim(xlim)
    }

    if(!is.null(highlight)){

        if(is(highlight, "list")){

            for(idx in seq_along(highlight)){
                tmp <- highlight[[idx]]
                p <- p + annotate("rect", 
                              xmin = tmp[1], xmax = tmp[2],
                              ymin = 0, ymax = Inf, 
                              fill = highlight_color, alpha = highlight_alpha)
            }

        }else{
            p <- p + annotate("rect", 
                              xmin = highlight[1], xmax = highlight[2],
                              ymin = 0, ymax = Inf, 
                              fill = highlight_color, alpha = highlight_alpha)
        }

    }

    if(add_coaccess){
        mat <- grange2mt(peak, weightCol)
        peak_cor_matrix <- cor(t(mat), method = "pearson")
        upper_tri <- upper.tri(peak_cor_matrix, diag = FALSE)

        if(is.null(coaccess_cor_threshold)){
            coaccess_cor_threshold <- 0.5
            cat("using coaccess_cor_threshold of", coaccess_cor_threshold, "\n")
        }

        high_cor_idx <- which(abs(peak_cor_matrix) > coaccess_cor_threshold & upper_tri, arr.ind = TRUE)
        coaccess_df <- data.frame(Peak1 = rownames(peak_cor_matrix)[high_cor_idx[, 1]],
                                 Peak2 = colnames(peak_cor_matrix)[high_cor_idx[, 2]],
                                 Cor = peak_cor_matrix[high_cor_idx])
        coaccess_df <- coaccess_df[order(coaccess_df$Cor,decreasing = TRUE),]
        if(is.null(coaccess_top_n)){
            coaccess_top_n <- 3
        }

        coaccess_df1 <- head(coaccess_df, n=coaccess_top_n)
        coaccess_df2 <- tail(coaccess_df, n=coaccess_top_n)
        coaccess_df <- rbind(coaccess_df1, coaccess_df2)
            

        peak1_data <- do.call(rbind, lapply(coaccess_df$Peak1, parse_peak))
        peak2_data <- do.call(rbind, lapply(coaccess_df$Peak2, parse_peak))
        df_processed <- data.frame(peak1_chr = peak1_data$chr,
                                   peak1_start = as.numeric(peak1_data$start),
                                   peak1_end = as.numeric(peak1_data$end),
                                   peak2_chr = peak2_data$chr,
                                   peak2_start = as.numeric(peak2_data$start),
                                   peak2_end = as.numeric(peak2_data$end),
                                   score = coaccess_df$Cor)
        df_processed$peak1_mid <- (df_processed$peak1_start + df_processed$peak1_end) / 2
        df_processed$peak2_mid <- (df_processed$peak2_start + df_processed$peak2_end) / 2
        df_processed$distance <- df_processed$peak2_start - df_processed$peak1_end
        
        peak1_end <- peak2_start <- NULL

        access_p <- ggplot(df_processed) +
                        geom_curve(aes(x = peak1_end, y = 0, 
                                    xend = peak2_start, yend = 0,color = score),
                                   curvature = curvature) +
                        scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                                            midpoint = 0,limits = c(-1, 1)) + 
                        coord_cartesian(xlim = xlim, ylim = c(-0.1, 0)) +
                        theme_classic() +
                        theme(axis.title = element_blank(),
                              axis.text = element_blank(),
                              axis.ticks = element_blank(),
                              panel.grid = element_blank(),
                              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
                              legend.position = coaccess_legend_pos,
                              legend.text = element_text(size = coaccess_legend_text_size),  
                              legend.title = element_text(size = coaccess_legend_title_size))
       
    }

    if(add_cluster_tree && add_coaccess){
        
        if(is.null(design)){
            design <- "AAAB
                       CCC#"
        }
        
        
        all_p <- plot_list(p, tree_p, access_p, design = design) &  theme(plot.margin=margin(0,0,0,0))

    }else if(add_cluster_tree){
        if(is.null(design)){
            design <- "AAAB"
        }

        all_p <- plot_list(p, tree_p, design = design) &  theme(plot.margin=margin(0,0,0,0))

    }else if(add_coaccess){

        if(is.null(design)){
            design <- "A
                       B"
        }

        
        all_p <- plot_list(p, access_p, design = design) &  theme(plot.margin=margin(0,0,0,0))

    }else{
        if(interactive){
            p <- girafe(ggobj = p, width_svg = width_svg,
                        height_svg = height_svg,
                        options = list(opts_hover(css = "fill: orange; stroke: black; stroke-width: 2px;"),
                                       opts_hover_inv(css = "opacity: 0.5;"),
                                       opts_zoom(min = .5, max = 5),
                                       opts_tooltip(css = "background-color: white; border: 1px solid black; padding: 5px; border-radius: 3px;")))

        }
        all_p <- p
    }
    
    return(all_p)
}

#' @import S4Vectors IRanges
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom magrittr %>%
getChrCov <- function(peak, weightCol, chrs, xlim, lower=1) {
    if (is(peak, "GRanges")) {
        peak.gr <- peak
    } else if (file.exists(peak)) {
        peak.gr <- readPeakFile(peak, as="GRanges")
    } else {
        stop("peak should be a GRanges object or a peak file...")
    }

    if ( is.null(weightCol)) {
        peak.cov <- coverage(peak.gr)
    } else {
        weight <- mcols(peak.gr)[[weightCol]]
        peak.cov <- coverage(peak.gr, weight=weight)
    }

    cov <- lapply(peak.cov, IRanges::slice, lower=lower)

    get.runValue <- function(x) {
        y <- runValue(x)
        sapply(y@listData, mean)
        ## value <- x@subject@values
        ## value[value != 0]
    }

    chr <- start <- end <- cnt <- NULL
    
    ldf <- lapply(1:length(cov), function(i) {
        x <- cov[[i]]
        if (length(x@ranges) == 0) {
            msg <- paste0(names(cov[i]),
                          " dosen't contain signal higher than ",
                          lower)
            message(msg)
            return(NA)
        }
        data.frame(chr   = names(cov[i]),
                   start = start(x),
                   end   = end(x),
                   cnt   = get.runValue(x)
                                        # the following versions are more slower
                                        # unlist(runValue(x)) 
                                        # sapply(x, runValue)
                   )
    })

    ldf <- ldf[!is.na(ldf)]
    df <- do.call("rbind", ldf)
    
    chr.sorted <- sortChrName(as.character(unique(df$chr)))
    df$chr <- factor(df$chr, levels=chr.sorted)
    if (!is.null(chrs) && !all(is.na(chrs)) && all(chrs %in% chr.sorted)) {
        df <- df[df$chr %in% chrs, ]
    }
    if (!is.null(xlim) && !all(is.na(xlim)) && is.numeric(xlim) && length(xlim) == 2) {
        df <- df[df$start >= xlim[1] & df$end <= xlim[2],]
    }

    df2 <- group_by(df, chr, start, end) %>% summarise(value=sum(cnt), .groups = "drop")
    return(df2)
}

# a simple `stringr::str_sort(numeric=TRUE)` implementation
sortChrName <- function(chr.name, decreasing = FALSE) {
    ## universal sort function, support organisms other than human
    chr_part <- sub("^(\\D*)(\\d*)$", "\\1", chr.name)
    num_part <- as.numeric(sub("^(\\D*)(\\d*)$", "\\2", chr.name))
    chr.name[order(chr_part, num_part, decreasing = decreasing)]
}


