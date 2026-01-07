#' @title plotCov
#' @details Plot peak coverage
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
#' @param cluster_dist_method method for calculate cluster tree. Details see [stats::dist()]
#' @param cluster_hclust_methond method for hclust. Details see [stats::hclust()]
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
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats as.formula
#' @importFrom S4Vectors head
#' @importFrom S4Vectors tail
#' @export
#' @examples
#' peakfile <- system.file("extdata", "sample_peaks.txt", package = "epiSeeker")
#' peak <- readPeakFile(peakfile)
#' plotCov(peak)
#' @author G Yu
plotCov <- function(peak, weightCol = NULL,
                    facet_level = NULL,
                    highlight = NULL,
                    highlight_color = "#c6c3c3",
                    highlight_alpha = 0.2,
                    xlab = "Chromosome Size (bp)",
                    ylab = "",
                    interactive = FALSE,
                    width_svg = 10,
                    height_svg = 6,
                    title = "ChIP Peaks over Chromosomes",
                    x_text_size = 10,
                    y_text_size = 10,
                    facet_label_text_size = 10,
                    chrs = NULL,
                    xlim = NULL,
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
                    coaccess_legend_pos = c(0.9, 0.5),
                    coaccess_legend_text_size = 10,
                    coaccess_legend_title_size = 12) {
    facet_var <- match.arg(facet_var, c("chr~.", ".~chr", ".~.id", ".id~.", ".id~chr", "chr~.id"))

    # 1. Get coverage data
    tm <- .get_cov_data(peak, weightCol, chrs, xlim, lower, facet_level)

    # 2. Add cluster tree if requested
    tree_res <- if (add_cluster_tree) {
        .handle_cluster_tree(peak, weightCol, cluster_dist_method, cluster_hclust_methond)
    } else {
        NULL
    }
    if (!is.null(tree_res)) {
        tm$.id <- factor(tm$.id, levels = tree_res$facet_level, labels = tree_res$facet_level)
    }

    # 3. Create main coverage plot
    p <- .create_cov_plot_base(tm, peak, fill_color, interactive, xlab, ylab, title, facet_var, facet_scales, facet_label_text_size, x_text_size, y_text_size, legend_position, xlim, highlight, highlight_color, highlight_alpha)

    # 4. Add co-accessibility or interactive wrapping
    all_p <- .finalize_all_plot(p, peak, weightCol, add_coaccess, add_cluster_tree, tree_res$tree_p, coaccess_cor_threshold, coaccess_top_n, curvature, xlim, coaccess_legend_pos, coaccess_legend_text_size, coaccess_legend_title_size, design, interactive, width_svg, height_svg)

    return(all_p)
}

.get_cov_data <- function(peak, weightCol, chrs, xlim, lower, facet_level) {
    isList <- is.list(peak)
    if (!isList) {
        tm <- getChrCov(peak = peak, weightCol = weightCol, chrs = chrs, xlim = xlim, lower = lower)
    } else {
        ltm <- lapply(peak, getChrCov, weightCol = weightCol, chrs = chrs, xlim = xlim, lower = lower)
        if (is.null(names(ltm))) {
            nn <- paste0("peak", seq_along(ltm))
            warning("input is not a named list, set the name automatically to ", paste(nn, collapse = " "))
            names(ltm) <- nn
        }
        tm <- dplyr::bind_rows(ltm, .id = ".id")
        chr.sorted <- sortChrName(as.character(unique(tm$chr)))
        tm$chr <- factor(tm$chr, levels = chr.sorted)
        if (!is.null(facet_level)) {
            tm$.id <- factor(tm$.id, levels = facet_level, labels = facet_level)
        }
    }
    return(tm)
}

.handle_cluster_tree <- function(peak, weightCol, cluster_dist_method, cluster_hclust_methond) {
    if (!is.list(peak)) stop("add_cluster_tree function apply to a list of peak")
    mat <- grange2mt(peak, weightCol)
    d <- dist(t(mat), method = cluster_dist_method)
    hc <- hclust(d, method = cluster_hclust_methond)
    rlang::check_installed(c("ape", "ggtree"))
    tree <- if (requireNamespace("ape", quietly = TRUE)) ape::as.phylo(hc) else NULL
    tree_p <- if (requireNamespace("ggtree", quietly = TRUE)) ggtree::ggtree(tree) + ggplot2::scale_x_reverse() else NULL
    tree_p_data <- as.data.frame(tree_p$data)
    tree_p_data <- tree_p_data[tree_p_data$isTip, ]
    facet_level <- tree_p_data[order(tree_p_data$y, decreasing = TRUE), "label"]
    list(tree_p = tree_p, facet_level = facet_level)
}

.create_cov_plot_base <- function(tm, peak, fill_color, interactive, xlab, ylab, title, facet_var, facet_scales, facet_label_text_size, x_text_size, y_text_size, legend_position, xlim, highlight, highlight_color, highlight_alpha) {
    if (length(tm$chr) == 0) {
        p <- ggplot(data.frame(x = 1)) + geom_blank()
        return(p)
    }
    chr <- start <- end <- value <- .id <- NULL
    p <- ggplot(tm, aes(.data$start, .data$value))
    isList <- is.list(peak)
    if (isList) {
        cols <- if (length(fill_color) == length(peak) && all(is_valid_color(fill_color))) fill_color else generate_colors(fill_color, n = length(peak))
        if (interactive && requireNamespace("ggiraph", quietly = TRUE)) {
            p <- p + ggiraph::geom_rect_interactive(aes(xmin = .data$start, ymin = 0, xmax = .data$end, ymax = .data$value, fill = .data$.id, color = .data$.id, tooltip = paste0("peak start: ", .data$start, "\n", "peak end: ", .data$end, "\n", "peak value: ", round(.data$value, 2)), data_id = paste0(chr, "_", .data$start, "_", .data$end)), hover_nearest = TRUE)
        } else {
            p <- p + geom_rect(aes(xmin = .data$start, ymin = 0, xmax = .data$end, ymax = .data$value, fill = .data$.id, color = .data$.id))
        }
        p <- p + scale_color_manual(values = cols) + scale_fill_manual(values = cols)
    } else {
        if (interactive && requireNamespace("ggiraph", quietly = TRUE)) {
            p <- p + ggiraph::geom_rect_interactive(aes(xmin = .data$start, ymin = 0, xmax = .data$end, ymax = .data$value, tooltip = paste0("peak start: ", .data$start, "\n", "peak end: ", .data$end, "\n", "peak value: ", round(.data$value, 2)), data_id = paste0(chr, "_", .data$start, "_", .data$end)), fill = fill_color, color = fill_color, hover_nearest = TRUE)
        } else {
            p <- p + geom_rect(aes(xmin = .data$start, ymin = 0, xmax = .data$end, ymax = .data$value), fill = fill_color, color = fill_color)
        }
    }

    if (length(unique(tm$chr)) > 1 || length(peak) > 1) {
        if (is.null(facet_var)) {
            facet_var <- if (length(unique(tm$chr)) > 1 && length(peak) > 1) "chr ~ .id" else if (length(unique(tm$chr)) > 1) "chr ~." else ".id ~."
        }
        p <- p + facet_grid(as.formula(facet_var), scales = facet_scales)
    }
    p <- p + theme_classic() + labs(x = xlab, y = ylab, title = title, fill = NULL, color = NULL) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(strip.text.y = element_text(angle = 360, size = facet_label_text_size), axis.text.x = element_text(size = x_text_size), axis.text.y = element_text(size = y_text_size)) +
        scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_si("")))

    if (!is.null(legend_position)) p <- p + theme(legend.position = legend_position)
    if (!is.null(xlim) && !all(is.na(xlim)) && is.numeric(xlim) && length(xlim) == 2) p <- p + xlim(xlim)
    if (!is.null(highlight)) {
        highlights <- if (is.list(highlight)) highlight else list(highlight)
        for (tmp in highlights) {
            p <- p + annotate("rect", xmin = tmp[1], xmax = tmp[2], ymin = 0, ymax = Inf, fill = highlight_color, alpha = highlight_alpha)
        }
    }
    return(p)
}

.finalize_all_plot <- function(p, peak, weightCol, add_coaccess, add_cluster_tree, tree_p, coaccess_cor_threshold, coaccess_top_n, curvature, xlim, coaccess_legend_pos, coaccess_legend_text_size, coaccess_legend_title_size, design, interactive, width_svg, height_svg) {
    access_p <- if (add_coaccess) .create_coaccess_plot(peak, weightCol, coaccess_cor_threshold, coaccess_top_n, curvature, xlim, coaccess_legend_pos, coaccess_legend_text_size, coaccess_legend_title_size) else NULL

    if (add_cluster_tree && add_coaccess) {
        design <- design %||% "AAAB\nCCC#"
        all_p <- plot_list(p, tree_p, access_p, design = design) & theme(plot.margin = margin(0, 0, 0, 0))
    } else if (add_cluster_tree) {
        design <- design %||% "AAAB"
        all_p <- plot_list(p, tree_p, design = design) & theme(plot.margin = margin(0, 0, 0, 0))
    } else if (add_coaccess) {
        design <- design %||% "A\nB"
        all_p <- plot_list(p, access_p, design = design) & theme(plot.margin = margin(0, 0, 0, 0))
    } else {
        if (interactive && requireNamespace("ggiraph", quietly = TRUE)) {
            p <- ggiraph::girafe(ggobj = p, width_svg = width_svg, height_svg = height_svg, options = list(ggiraph::opts_hover(css = "fill: orange; stroke: black; stroke-width: 2px;"), ggiraph::opts_hover_inv(css = "opacity: 0.5;"), ggiraph::opts_zoom(min = .5, max = 5), ggiraph::opts_tooltip(css = "background-color: white; border: 1px solid black; padding: 5px; border-radius: 3px;")))
        }
        all_p <- p
    }
    return(all_p)
}

.create_coaccess_plot <- function(peak, weightCol, coaccess_cor_threshold, coaccess_top_n, curvature, xlim, coaccess_legend_pos, coaccess_legend_text_size, coaccess_legend_title_size) {
    mat <- grange2mt(peak, weightCol)
    peak_cor_matrix <- cor(t(mat), method = "pearson")
    upper_tri <- upper.tri(peak_cor_matrix, diag = FALSE)
    coaccess_cor_threshold <- coaccess_cor_threshold %||% 0.5
    high_cor_idx <- which(abs(peak_cor_matrix) > coaccess_cor_threshold & upper_tri, arr.ind = TRUE)
    coaccess_df <- data.frame(Peak1 = rownames(peak_cor_matrix)[high_cor_idx[, 1]], Peak2 = colnames(peak_cor_matrix)[high_cor_idx[, 2]], Cor = peak_cor_matrix[high_cor_idx])
    coaccess_df <- coaccess_df[order(coaccess_df$Cor, decreasing = TRUE), ]
    coaccess_top_n <- coaccess_top_n %||% 3
    coaccess_df <- rbind(head(coaccess_df, n = coaccess_top_n), tail(coaccess_df, n = coaccess_top_n))

    peak1_data <- do.call(rbind, lapply(coaccess_df$Peak1, parse_peak))
    peak2_data <- do.call(rbind, lapply(coaccess_df$Peak2, parse_peak))
    df_processed <- data.frame(peak1_end = as.numeric(peak1_data$end), peak2_start = as.numeric(peak2_data$start), score = coaccess_df$Cor)
    score <- peak1_end <- peak2_start <- NULL
    ggplot(df_processed) +
        geom_curve(aes(x = .data$peak1_end, y = 0, xend = .data$peak2_start, yend = 0, color = .data$score), curvature = curvature) +
        scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
        coord_cartesian(xlim = xlim, ylim = c(-0.1, 0)) +
        theme_classic() +
        theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = coaccess_legend_pos, legend.text = element_text(size = coaccess_legend_text_size), legend.title = element_text(size = coaccess_legend_title_size))
}

#' @import IRanges
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom magrittr %>%
#' @importFrom S4Vectors runValue
getChrCov <- function(peak, weightCol, chrs, xlim, lower = 1) {
    if (is(peak, "GRanges")) {
        peak.gr <- peak
    } else if (file.exists(peak)) {
        peak.gr <- readPeakFile(peak, as = "GRanges")
    } else {
        stop("peak should be a GRanges object or a peak file...")
    }

    if (is.null(weightCol)) {
        peak.cov <- coverage(peak.gr)
    } else {
        weight <- mcols(peak.gr)[[weightCol]]
        peak.cov <- coverage(peak.gr, weight = weight)
    }

    cov <- lapply(peak.cov, IRanges::slice, lower = lower)

    get.runValue <- function(x) {
        y <- runValue(x)
        # sapply(y@listData, mean)
        vapply(y@listData, mean, FUN.VALUE = numeric(1))
        ## value <- x@subject@values
        ## value[value != 0]
    }

    chr <- start <- end <- cnt <- NULL

    ldf <- lapply(seq_len(length(cov)), function(i) {
        x <- cov[[i]]
        if (length(x@ranges) == 0) {
            msg <- paste0(
                names(cov[i]),
                " dosen't contain signal higher than ",
                lower
            )
            message(msg)
            return(NA)
        }
        data.frame(
            chr = names(cov[i]),
            start = start(x),
            end = end(x),
            cnt = get.runValue(x)
            # the following versions are more slower
            # unlist(runValue(x))
            # sapply(x, runValue)
        )
    })

    ldf <- ldf[!is.na(ldf)]
    df <- do.call("rbind", ldf)

    chr.sorted <- sortChrName(as.character(unique(df$chr)))
    df$chr <- factor(df$chr, levels = chr.sorted)
    if (!is.null(chrs) && !all(is.na(chrs)) && all(chrs %in% chr.sorted)) {
        df <- df[df$chr %in% chrs, ]
    }
    if (!is.null(xlim) && !all(is.na(xlim)) && is.numeric(xlim) && length(xlim) == 2) {
        df <- df[df$start >= xlim[1] & df$end <= xlim[2], ]
    }

    df2 <- group_by(df, chr, start, end) %>% summarise(value = sum(cnt), .groups = "drop")
    return(df2)
}

# a simple `stringr::str_sort(numeric=TRUE)` implementation
sortChrName <- function(chr.name, decreasing = FALSE) {
    ## universal sort function, support organisms other than human
    chr_part <- sub("^(\\D*)(\\d*)$", "\\1", chr.name)
    num_part <- as.numeric(sub("^(\\D*)(\\d*)$", "\\2", chr.name))
    chr.name[order(chr_part, num_part, decreasing = decreasing)]
}
