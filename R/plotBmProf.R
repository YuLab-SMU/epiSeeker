#' Plot base modification profile
#'
#' @title plotBmProf
#' @param df the base modification data.frame
#' @param motif_color the color for different motifs(CHH,CHG,CG)
#' @param title the title of the plot, can also be a list of title.
#' @param interactive produce interactive fig or not.
#' @param width_svg  width_svg.
#' @param height_svg height_svg.
#' @param xlim the specified interval of region, must be the sub-interval of the dmR. list for list df
#' @param highlight a region or a list of region to highlight.
#' @param highlight_color colors of highlight rect. Default "#c6c3c3"
#' @param highlight_alpha alpha of highlight rect.
#' @param xlab the x label, can also be a list of x label
#' @param ylab the y label, can also be a list of y label
#' @param second_ylab the ylab for second y-axis
#' @param switch_y_value switch the value from left y-axis to right y-axis
#' @param legend_lab_motif the label of legend for motif
#' @param legend_lab_value2 the label of legend for the second value(ylab is the label for the first value)
#' @param strip_placement strip.placement
#' @param angle_of_facet_label the angle of facet label, e.g. 0 is horizontal
#' @param alpha transparency for the depth information line
#' @param y_ticks_length the length of y-axis ticks
#' @param x_ticks_length the length of x-axis ticks
#' @param auto_x_axis use auto x axis or not.
#' @param strip_border add border to the facet label or not
#' @param facet_label_text_size the size of facet label text
#' @param axis_title_text_size the size of axis title text
#' @param title_text_size the size of the title text
#' @param right_y_axis_text_size the size of the left y axis text,this work when depth information is taken into account
#' @param left_y_axis_text_size the size of the left y axis text
#' @param x_axis_text_size the size of x axis text
#' @param depth_heatmap draw the heatmap of depth information or not
#' @param nrow the nrow of plotting a list of dmR
#' @param ncol the ncol of plotting a list of dmR
#' @param panel_spacing the distance between panels
#' @param legend_box_spacing the distance between legend and plotting area,"cm"
#' @param legend_position the position of legend
#' @return ggplot object
#' @importFrom aplot plot_list
#' @importFrom methods is
#' @examples
#' require(BSgenome.Hsapiens.UCSC.hg38)
#' data(demo_bmdata)
#' bmMatrix <- getBmMatrix(
#'     region = data.frame(chr = "chr22", start = 10525991, end = 10526342),
#'     BSgenome = BSgenome.Hsapiens.UCSC.hg38,
#'     input = demo_bmdata,
#'     #                          base = "C",
#'     motif = c("CG")
#' )
#' plotBmProf(bmMatrix)
#' @export
plotBmProf <- function(df,
                       motif_color = NULL,
                       title = NULL,
                       xlim = NULL,
                       interactive = FALSE,
                       width_svg = 10,
                       height_svg = 6,
                       highlight = NULL,
                       highlight_color = "#c6c3c3",
                       highlight_alpha = 0.2,
                       xlab = "Genomic Region(5'->3')",
                       ylab = NULL,
                       second_ylab = NULL,
                       switch_y_value = TRUE,
                       legend_lab_motif = NULL,
                       legend_lab_value2 = NULL,
                       strip_placement = "outside",
                       angle_of_facet_label = 360,
                       alpha = 0.6,
                       y_ticks_length = 0.25,
                       x_ticks_length = 0.25,
                       auto_x_axis = TRUE,
                       strip_border = FALSE,
                       facet_label_text_size = 10,
                       axis_title_text_size = 17,
                       title_text_size = 20,
                       right_y_axis_text_size = 10,
                       left_y_axis_text_size = 10,
                       x_axis_text_size = 10,
                       depth_heatmap = TRUE,
                       nrow = NULL,
                       ncol = NULL,
                       panel_spacing = 1,
                       legend_box_spacing = 3,
                       legend_position = "right") {
    # assign default value for label
    if (is(df, "list")) {
        tmpdf <- df[[1]]
    } else {
        tmpdf <- df
    }

    vName <- unique(tmpdf$type)
    n0 <- length(vName)


    if (is.null(legend_lab_motif)) {
        legend_lab_motif <- "motif"
    }

    if (n0 == 1) {
        if (is.null(ylab)) ylab <- vName
    } else {
        vName1 <- vName[1]
        vName2 <- vName[2]

        if (switch_y_value) {
            vName1 <- vName[2]
            vName2 <- vName[1]
        }

        if (is.null(ylab)) ylab <- vName1
        if (is.null(second_ylab)) second_ylab <- vName2
        if (is.null(legend_lab_value2)) legend_lab_value2 <- vName2
    }

    if (is.null(nrow) && is.null(ncol)) {
        ncol <- 1
    }

    if (is(df, "list")) {
        # check xlab,ylab and title
        if (length(xlab) != 1) {
            if (length(xlab) != length(df)) {
                stop("please input xlab with correct length...")
            }
        } else {
            xlab <- rep(xlab, length(df))
        }

        if (length(ylab) != 1) {
            if (length(ylab) != length(df)) {
                stop("please input ylab with correct length...")
            }
        } else {
            ylab <- rep(ylab, length(df))
        }


        # check xlim
        if (!is.null(xlim)) {
            if (length(xlim != length(df))) {
                stop("the length of xlim and the length of df are not equal...")
            }
        }

        if (!is.null(title)) {
            if (length(title) != length(df)) {
                stop("please input title with correct length...")
            }

            temp <- list()

            for (i in seq_len(length(df))) {
                p <- plotBmProf.internal(
                    df = df[[i]],
                    motif_color = motif_color,
                    title = title[i],
                    xlim = xlim[[i]],
                    interactive = interactive,
                    width_svg = width_svg,
                    height_svg = height_svg,
                    highlight = highlight,
                    highlight_color = highlight_color,
                    highlight_alpha = highlight_alpha,
                    xlab = xlab[i],
                    ylab = ylab[i],
                    second_ylab = second_ylab,
                    switch_y_value = switch_y_value,
                    legend_lab_motif = legend_lab_motif,
                    legend_lab_value2 = legend_lab_value2,
                    strip_placement = strip_placement,
                    angle_of_facet_label = angle_of_facet_label,
                    alpha = alpha,
                    y_ticks_length = y_ticks_length,
                    x_ticks_length = x_ticks_length,
                    auto_x_axis = auto_x_axis,
                    strip_border = strip_border,
                    facet_label_text_size = facet_label_text_size,
                    axis_title_text_size = axis_title_text_size,
                    title_text_size = title_text_size,
                    right_y_axis_text_size = right_y_axis_text_size,
                    left_y_axis_text_size = left_y_axis_text_size,
                    x_axis_text_size = x_axis_text_size,
                    depth_heatmap = depth_heatmap,
                    panel_spacing = panel_spacing,
                    legend_box_spacing = legend_box_spacing,
                    legend_position = legend_position
                )

                temp[[i]] <- p
            }

            p <- plot_list(
                gglist = temp,
                ncol = ncol,
                nrow = nrow
            )
            return(p)
        }

        temp <- list()

        for (i in seq_len(length(df))) {
            p <- plotBmProf.internal(
                df = df[[i]],
                motif_color = motif_color,
                title = title,
                xlim = xlim[[i]],
                highlight = highlight[[i]],
                highlight_color = highlight_color,
                highlight_alpha = highlight_alpha,
                interactive = interactive,
                width_svg = width_svg,
                height_svg = height_svg,
                xlab = xlab[i],
                ylab = ylab[i],
                second_ylab = second_ylab,
                switch_y_value = switch_y_value,
                legend_lab_motif = legend_lab_motif,
                legend_lab_value2 = legend_lab_value2,
                strip_placement = strip_placement,
                angle_of_facet_label = angle_of_facet_label,
                alpha = alpha,
                y_ticks_length = y_ticks_length,
                x_ticks_length = x_ticks_length,
                auto_x_axis = auto_x_axis,
                strip_border = strip_border,
                facet_label_text_size = facet_label_text_size,
                axis_title_text_size = axis_title_text_size,
                title_text_size = title_text_size,
                right_y_axis_text_size = right_y_axis_text_size,
                left_y_axis_text_size = left_y_axis_text_size,
                x_axis_text_size = x_axis_text_size,
                depth_heatmap = depth_heatmap,
                panel_spacing = panel_spacing,
                legend_box_spacing = legend_box_spacing,
                legend_position = legend_position
            )

            temp[[i]] <- p
        }


        p <- plot_list(
            gglist = temp,
            ncol = ncol,
            nrow = nrow
        )
    } else {
        p <- plotBmProf.internal(
            df = df,
            motif_color = motif_color,
            title = title,
            xlim = xlim,
            highlight = highlight,
            highlight_color = highlight_color,
            highlight_alpha = highlight_alpha,
            xlab = xlab,
            ylab = ylab,
            interactive = interactive,
            width_svg = width_svg,
            height_svg = height_svg,
            second_ylab = second_ylab,
            switch_y_value = switch_y_value,
            legend_lab_motif = legend_lab_motif,
            legend_lab_value2 = legend_lab_value2,
            strip_placement = strip_placement,
            angle_of_facet_label = angle_of_facet_label,
            alpha = alpha,
            y_ticks_length = y_ticks_length,
            x_ticks_length = x_ticks_length,
            auto_x_axis = auto_x_axis,
            strip_border = strip_border,
            facet_label_text_size = facet_label_text_size,
            axis_title_text_size = axis_title_text_size,
            title_text_size = title_text_size,
            right_y_axis_text_size = right_y_axis_text_size,
            left_y_axis_text_size = left_y_axis_text_size,
            x_axis_text_size = x_axis_text_size,
            depth_heatmap = depth_heatmap,
            panel_spacing = panel_spacing,
            legend_box_spacing = legend_box_spacing,
            legend_position = legend_position
        )
    }


    return(p)
}


#' @import ggplot2
#' @importFrom scales rescale
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges seqnames
#' @importFrom IRanges IRanges
#' @importFrom aplot xlim2
#' @importFrom aplot insert_bottom
#' @importFrom magrittr %>%
plotBmProf.internal <- function(df,
                                motif_color,
                                interactive = FALSE,
                                width_svg = 10,
                                height_svg = 6,
                                title,
                                xlim,
                                highlight,
                                highlight_color,
                                highlight_alpha,
                                xlab,
                                ylab,
                                second_ylab,
                                switch_y_value,
                                legend_lab_motif,
                                legend_lab_value2,
                                strip_placement = "outside",
                                angle_of_facet_label = 0,
                                alpha = 0.3,
                                y_ticks_length = 0.25,
                                x_ticks_length = 0.25,
                                auto_x_axis = TRUE,
                                strip_border = FALSE,
                                facet_label_text_size = 12,
                                axis_title_text_size = 17,
                                title_text_size = 20,
                                right_y_axis_text_size = 13,
                                left_y_axis_text_size = 13,
                                x_axis_text_size = 13,
                                depth_heatmap = TRUE,
                                panel_spacing = 1,
                                legend_box_spacing = 3,
                                legend_position = "right") {
    # 1. Prepare data and coordinate info
    prep <- .prepare_bm_data(df, xlim, switch_y_value)
    df <- prep$df
    coordinate <- prep$coordinate
    vName1 <- prep$vName1
    vName2 <- prep$vName2
    value1_max <- prep$value1_max
    value2_max <- prep$value2_max
    df_chr <- prep$df_chr
    n0 <- prep$n0

    # 2. Add main plot layers (methylation and depth)
    p <- .add_bm_layers(df, n0, vName1, vName2, value1_max, value2_max, interactive, legend_lab_motif, legend_lab_value2, second_ylab, alpha)

    # 3. Add highlights
    p <- .add_highlights(p, highlight, n0, value1_max, value2_max, highlight_color, highlight_alpha)

    # 4. Final themes and labels
    p <- .add_themes_and_labels(
        p, df, title, xlab, ylab, second_ylab, df_chr, coordinate,
        angle_of_facet_label, facet_label_text_size, strip_border,
        auto_x_axis, title_text_size, axis_title_text_size,
        left_y_axis_text_size, right_y_axis_text_size, x_axis_text_size,
        motif_color, y_ticks_length, x_ticks_length, panel_spacing,
        legend_box_spacing, legend_position, interactive, width_svg, height_svg
    )

    return(p)
}

.prepare_bm_data <- function(df, xlim, switch_y_value) {
    coordinate <- unique(df$coordinate)
    df_chr <- attr(df, "chromosome")
    vName <- unique(df$type)
    n0 <- length(vName)

    if (n0 == 1) {
        vName1 <- vName
        vName2 <- NULL
        value1_max <- ceiling(max(df$value))
        value2_max <- NULL
    } else {
        vName1 <- if (switch_y_value) vName[2] else vName[1]
        vName2 <- if (switch_y_value) vName[1] else vName[2]
        value1_max <- ceiling(max(df$value[df$type == vName1]))
        value2_max <- ceiling(max(df$value[df$type == vName2]))
    }

    df$value[df$type == vName1 & df$strand == "-"] <- (-1) * df$value[df$type == vName1 & df$strand == "-"]

    if (!is.null(xlim)) {
        coordinate <- c(xlim[1]:xlim[2])
        if (xlim[1] < min(df$coordinate) || xlim[2] > max(df$coordinate)) {
            df <- df[df$coordinate %in% coordinate, ]
        } else {
            new_df <- list()
            for (type_idx in unique(df$type)) {
                for (sample_idx in unique(df$sample)) {
                    tmp <- df[df$type == type_idx & df$sample == sample_idx, ]
                    tmp_df <- data.frame(coordinate = xlim[1]:xlim[2])
                    tmp_merge <- merge(tmp, tmp_df, all = TRUE)
                    tmp_merge$motif[is.na(tmp_merge$motif)] <- "none"
                    tmp_merge$strand[is.na(tmp_merge$strand)] <- "*"
                    tmp_merge$sample[is.na(tmp_merge$sample)] <- sample_idx
                    tmp_merge$value[is.na(tmp_merge$value)] <- 0
                    tmp_merge$type[is.na(tmp_merge$type)] <- type_idx
                    new_df[[paste0(type_idx, sample_idx)]] <- tmp_merge
                }
            }
            df <- do.call("rbind", new_df)
        }
    }

    # replace the "none" with unique(df$motif)[1]
    none_fill_in <- unique(df$motif[df$motif != "none"])[1]
    df$motif[df$motif == "none"] <- none_fill_in

    list(
        df = df, coordinate = coordinate, vName1 = vName1, vName2 = vName2,
        value1_max = value1_max, value2_max = value2_max, df_chr = df_chr, n0 = n0
    )
}

.add_bm_layers <- function(df, n0, vName1, vName2, value1_max, value2_max, interactive, legend_lab_motif, legend_lab_value2, second_ylab, alpha) {
    if (n0 == 2) {
        res <- .calculate_rescaled_values(df, vName2, value1_max, value2_max)
        pos_strand <- res$pos_strand
        neg_strand <- res$neg_strand

        if (interactive) {
            rlang::check_installed("ggiraph", reason = "For interactive plot.")
            p <- ggplot() +
                ggiraph::geom_col_interactive(
                    data = df[df$type == vName1, ],
                    mapping = aes(
                        x = .data$coordinate, y = .data$value, fill = .data$motif, color = .data$motif,
                        tooltip = paste0("Coordinate: ", .data$coordinate, "\nValue: ", round(.data$value, 2)),
                        data_id = .data$coordinate
                    ), size = 1
                )
        } else {
            p <- ggplot(df) +
                geom_col(data = df[df$type == vName1, ], mapping = aes(x = .data$coordinate, y = .data$value, fill = .data$motif, color = .data$motif))
        }

        p <- p +
            geom_line(data = pos_strand, mapping = aes(x = .data$coordinate, y = .data$value, linetype = paste0(vName2, " information")), color = "#868686FF", alpha = alpha) +
            geom_line(data = neg_strand, mapping = aes(x = .data$coordinate, y = .data$value, linetype = paste0(vName2, " information")), color = "#868686FF", alpha = alpha) +
            labs(fill = legend_lab_motif, linetype = legend_lab_value2) +
            guides(fill = guide_legend(order = 1), linetype = guide_legend(order = 0), color = "none")

        # reorganized axis
        if (nrow(pos_strand) == 1) {
            p <- p + scale_y_continuous(
                sec.axis = sec_axis(transform = ~ . * (value2_max / value1_max), name = second_ylab, breaks = c(-value2_max, 0), labels = c(value2_max, 0)),
                breaks = c(-value1_max, -0.5 * value1_max, 0), labels = c(paste0("(3'->5') ", value1_max), 0.5 * value1_max, 0)
            )
        } else if (nrow(neg_strand) == 1) {
            p <- p + scale_y_continuous(
                sec.axis = sec_axis(transform = ~ . * (value2_max / value1_max), name = second_ylab, breaks = c(0, value2_max), labels = c(0, value2_max)),
                breaks = c(-value1_max, -0.5 * value1_max, 0), labels = c(paste0("(3'->5') ", value1_max), 0.5 * value1_max, 0)
            )
        } else {
            p <- p + scale_y_continuous(
                sec.axis = sec_axis(transform = ~ rescale(., c(-value2_max, value2_max)), name = second_ylab, breaks = c(-value2_max, 0, value2_max), labels = c(value2_max, 0, value2_max)),
                breaks = c((-1) * value1_max, (-0.5) * value1_max, 0, (0.5) * value1_max, value1_max),
                labels = c(paste0("(3'->5') ", value1_max), (0.5) * value1_max, 0, (0.5) * value1_max, paste0("(5'->3') ", value1_max))
            )
        }
    } else {
        if (interactive) {
            rlang::check_installed("ggiraph", reason = "For interactive plot.")
            p <- ggplot(df) +
                ggiraph::geom_col_interactive(
                    mapping = aes(
                        x = .data$coordinate, y = .data$value, fill = .data$motif, color = .data$motif,
                        tooltip = paste0("Coordinate: ", .data$coordinate, "\nValue: ", round(.data$value, 2)),
                        data_id = .data$coordinate
                    ), size = 1
                )
        } else {
            p <- ggplot(df) +
                geom_col(mapping = aes(x = .data$coordinate, y = .data$value, fill = .data$motif, color = .data$motif))
        }
        p <- p + labs(fill = legend_lab_motif) + guides(color = "none") +
            scale_y_continuous(
                breaks = c((-1) * value1_max, (-0.5) * value1_max, 0, (0.5) * value1_max, value1_max),
                labels = c(paste0("(3'->5') ", value1_max), (0.5) * value1_max, 0, (0.5) * value1_max, paste0("(5'->3') ", value1_max))
            )
    }
    return(p)
}

.calculate_rescaled_values <- function(df, vName2, value1_max, value2_max) {
    pos_strand <- neg_strand <- df[df$type == vName2, ]
    pos_strand$value[pos_strand$strand == "-"] <- 0
    neg_strand$value[neg_strand$strand == "+"] <- 0

    ncol_tmp <- nrow(pos_strand)
    pos_val <- pos_strand$value
    neg_val <- neg_strand$value
    pos_val[ncol_tmp + 1] <- neg_val[ncol_tmp + 1] <- 0
    pos_val[ncol_tmp + 2] <- neg_val[ncol_tmp + 2] <- value2_max

    pos_strand$value <- rescale(pos_val, c(0, value1_max))[seq_len(ncol_tmp)]
    neg_strand$value <- rescale(neg_val, c(0, (-1) * value1_max))[seq_len(ncol_tmp)]

    pos_strand <- if (any(pos_strand$value != 0)) {
        pos_strand[pos_strand$value != 0, ]
    } else {
        tmp <- df[1, ]
        tmp$value <- 0
        tmp
    }
    neg_strand <- if (any(neg_strand$value != 0)) {
        neg_strand[neg_strand$value != 0, ]
    } else {
        tmp <- df[1, ]
        tmp$value <- 0
        tmp
    }

    list(pos_strand = pos_strand, neg_strand = neg_strand)
}

.add_highlights <- function(p, highlight, n0, value1_max, value2_max, highlight_color, highlight_alpha) {
    if (is.null(highlight)) {
        return(p)
    }
    y_max <- if (n0 == 1) value1_max else value2_max
    highlights <- if (is.list(highlight)) highlight else list(highlight)

    for (tmp in highlights) {
        p <- p + annotate("rect", xmin = tmp[1], xmax = tmp[2], ymin = -y_max, ymax = y_max, fill = highlight_color, alpha = highlight_alpha)
    }
    p <- gginnards::move_layers(x = p, idx = 4L, position = "bottom")
    return(p)
}

.add_themes_and_labels <- function(p, df, title, xlab, ylab, second_ylab, df_chr, coordinate,
                                   angle_of_facet_label, facet_label_text_size, strip_border,
                                   auto_x_axis, title_text_size, axis_title_text_size,
                                   left_y_axis_text_size, right_y_axis_text_size, x_axis_text_size,
                                   motif_color, y_ticks_length, x_ticks_length, panel_spacing,
                                   legend_box_spacing, legend_position, interactive, width_svg, height_svg) {
    p <- p + theme_classic()

    if (length(unique(df$sample)) > 1) {
        p <- p + facet_grid(sample ~ .) +
            theme(strip.text.y = element_text(angle = angle_of_facet_label, color = "black", face = "bold", size = facet_label_text_size))
    }

    if (!strip_border) {
        p <- p + theme(strip.background = element_blank())
    }

    p <- p + theme(axis.line.x = element_blank()) + geom_hline(yintercept = 0)

    if (!auto_x_axis) {
        len <- length(coordinate)
        breaks <- coordinate[c(1, floor(len * 0.25), floor(len * 0.5), floor(len * 0.75), len)]
        labels <- paste0(df_chr, ":", breaks)
        p <- p + scale_x_continuous(breaks = breaks, labels = labels)
    }

    title <- title %||% paste0(df_chr, ": ", coordinate[1], "~", coordinate[length(coordinate)], " bp")

    p <- p + labs(title = title, x = xlab, y = ylab) +
        theme(
            plot.title = element_text(hjust = 0.5, size = title_text_size),
            axis.title = element_text(size = axis_title_text_size),
            axis.text.x = element_text(vjust = -1.2, size = x_axis_text_size),
            axis.title.x = element_text(vjust = -1.2),
            axis.title.y.right = element_text(vjust = 5.5),
            axis.text.y.left = element_text(size = left_y_axis_text_size),
            axis.text.y.right = element_text(size = right_y_axis_text_size)
        )

    if (!is.null(motif_color)) {
        if (length(motif_color) != length(unique(df$motif))) stop("the length of color and the motif should be equal...")
        p <- p + scale_fill_manual(values = motif_color)
    }

    p <- p + theme(
        axis.ticks.length.y = unit(y_ticks_length, "cm"),
        axis.ticks.length.x = unit(x_ticks_length, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.spacing.y = unit(panel_spacing, "cm"),
        legend.box_spacing = unit(legend_box_spacing, "cm"),
        legend.position = legend_position
    )

    if (interactive && requireNamespace("ggiraph", quietly = TRUE)) {
        p <- ggiraph::girafe(
            ggobj = p, width_svg = width_svg, height_svg = height_svg,
            options = list(
                ggiraph::opts_hover(css = "fill: orange; stroke: black; stroke-width: 2px;"),
                ggiraph::opts_hover_inv(css = "opacity: 0.5;"),
                ggiraph::opts_zoom(min = .5, max = 5),
                ggiraph::opts_tooltip(css = "background-color: white; border: 1px solid black; padding: 5px; border-radius: 3px;")
            )
        )
    }
    return(p)
}
