#' plot base modification profile
#'
#'
#' @title plotBmProf
#' @param df the base modification dataframe
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
                       legend_position = "right"){

  # assign default value for label
  if(is(df,"list")){
    tmpdf <- df[[1]]
  }else{
    tmpdf <- df
  }

  vName <- unique(tmpdf$type)
  n0 <- length(vName)


  if(is.null(legend_lab_motif )){
    legend_lab_motif <- "motif"
  }

  if(n0 == 1){
    if(is.null(ylab)) ylab <- vName
  }else{

    vName1 <- vName[1]
    vName2 <- vName[2]

    if(switch_y_value){
      vName1 <- vName[2]
      vName2 <- vName[1]
    }

    if(is.null(ylab)) ylab <- vName1
    if(is.null(second_ylab)) second_ylab <- vName2
    if(is.null(legend_lab_value2)) legend_lab_value2 <- vName2
  }

  if(is.null(nrow) && is.null(ncol)){
    ncol <- 1
  }

  if(is(df,"list")){

    # check xlab,ylab and title
    if(length(xlab) != 1){

      if(length(xlab) != length(df)){
        stop("please input xlab with correct length...")
      }

    }else{

      xlab <- rep(xlab, length(df))

    }

    if(length(ylab) != 1){

      if(length(ylab) != length(df)){
        stop("please input ylab with correct length...")
      }

    }else{

      ylab <- rep(ylab, length(df))

    }


    # check xlim
    if(!is.null(xlim)){
      if(length(xlim != length(df))){
        stop("the length of xlim and the length of df are not equal...")
      }
    }

    if(!is.null(title)){

      if(length(title) != length(df)){
        stop("please input title with correct length...")
      }

      temp <- list()

      for (i in seq_len(length(df))) {

        p <- plotBmProf.internal(df = df[[i]],
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
                                 legend_position = legend_position)

        temp[[i]] <- p
      }

      p <- plot_list(gglist = temp,
                     ncol = ncol,
                     nrow = nrow)
      return(p)

    }

    temp <- list()

    for (i in seq_len(length(df))) {

      p <- plotBmProf.internal(df = df[[i]],
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
                               legend_position = legend_position)

      temp[[i]] <- p
    }


    p <- plot_list(gglist = temp,
                   ncol = ncol,
                   nrow = nrow)


  }else{
    p <- plotBmProf.internal(df = df,
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
                             legend_position = legend_position)
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
#' @importFrom gginnards move_layers
#' @importFrom ggiraph geom_col_interactive
#' @importFrom ggiraph geom_line_interactive
#' @importFrom ggiraph girafe
#' @importFrom ggiraph opts_hover
#' @importFrom ggiraph opts_hover_inv
#' @importFrom ggiraph opts_tooltip
#' @importFrom ggiraph opts_zoom
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
                                legend_position = "right"){

  # extract coordinate information
  coordinate <- unique(df$coordinate)
  df_chr <- attr(df,"chromosome")

  vName <- unique(df$type)
  n0 <- length(vName)

  if(n0 == 1){
    vName1 <- vName
    value1_max <- ceiling(max(df$value))

  }else{
    vName1 <- vName[1]
    vName2 <- vName[2]

    if(switch_y_value){
      vName1 <- vName[2]
      vName2 <- vName[1]
    }

    value1_max <- ceiling(max(df$value[df$type == vName1]))
    value2_max <- ceiling(max(df$value[df$type == vName2]))
  }

  # flip the value on minus strand
  df$value[df$type == vName1 & df$strand == "-"] <- (-1)*df$value[df$type == vName1 & df$strand == "-"]


  if(!is.null(xlim)){

    coordinate <- c(xlim[1]:xlim[2])
    # xlim must be a subinterval of the dmR
    if(xlim[1]<min(coordinate) || xlim[2]>max(coordinate)){
      df <- df[df$coordinate %in% coordinate,]
    }else{
      new_df <- list()
      for(type_idx in unique(df$type)){
          for(sample_idx in unique(df$sample)){
              tmp <- df[df$type == type_idx & df$sample == sample_idx,]
              tmp_df <- data.frame(coordinate = xlim[1]:xlim[2])
              tmp_merge <- merge(tmp, tmp_df, all = TRUE)
              tmp_merge$motif[is.na(tmp_merge$motif)] <- 'none'
              tmp_merge$strand[is.na(tmp_merge$strand)] <- '*'
              tmp_merge$sample[is.na(tmp_merge$sample)] <- sample_idx
              tmp_merge$value[is.na(tmp_merge$value)] <- 0
              tmp_merge$type[is.na(tmp_merge$type)] <- type_idx

              new_df[[paste0(type_idx, sample_idx)]] <- tmp_merge
          }
      }
      df <- do.call("rbind", new_df)
    }
    
    cat(">> plotting region from ",paste0(df_chr,":",xlim[1]),
        " to ",paste0(df_chr,":",xlim[2]),"...\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n",sep = "")

  }else{
    cat(">> plotting region from ",paste0(df_chr,":",coordinate[1]),
        " to ",paste0(df_chr,":",coordinate[length(coordinate)]),"...\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n",sep = "")
  }


  # replace the "none" with unique(df$motif)[1] for the convenience of plotting
  # every information in "none" item is 0, so there is no effect
  none_fill_in <- unique(df$motif[df$motif != "none"])[1]
  df$motif[df$motif == "none"] <- none_fill_in

  # extract the tiltle information
  if(is.null(title)){
    title <- paste0(df_chr,
                    ": ",
                    coordinate[1],
                    "~",
                    coordinate[length(coordinate)],
                    " bp")
  }


  # global binding for value and motif
  value <- motif <- NULL

  if(n0 == 2){

    positive_strand_temp <- negative_strand_temp <- df[df$type==vName2,]
    positive_strand_temp$value[positive_strand_temp$strand == "-"] <- 0
    negative_strand_temp$value[negative_strand_temp$strand == "+"] <- 0

    ncol_tmp <- nrow(positive_strand_temp)
    positive_strand_temp_value <- positive_strand_temp$value
    negative_strand_temp_value <- negative_strand_temp$value

    # add value to correct the rescale
    positive_strand_temp_value[ncol_tmp+1] <- negative_strand_temp_value[ncol_tmp+1] <- 0
    positive_strand_temp_value[ncol_tmp+2] <- negative_strand_temp_value[ncol_tmp+2] <- value2_max

    rescale_positive_strand <- rescale(positive_strand_temp_value,c(0,value1_max))[1:ncol_tmp]
    rescale_negative_strand <- rescale(negative_strand_temp_value,c(0,(-1)*value1_max))[1:ncol_tmp]

    positive_strand_temp$value <- rescale_positive_strand
    negative_strand_temp$value <- rescale_negative_strand

    positive_strand_temp <- positive_strand_temp[positive_strand_temp$value != 0,]
    negative_strand_temp <- negative_strand_temp[negative_strand_temp$value != 0,]

    if(nrow(positive_strand_temp) == 0){
      positive_strand_temp <- df[1,]
      positive_strand_temp$value <- 0
    }

    if(nrow(negative_strand_temp) == 0){
      negative_strand_temp$coordinate <- df[1,]
      negative_strand_temp$value <- 0
    }

    # plot the methylation information and cover depth information
    if(interactive){

      # p <- ggplot() +
      #       geom_col_interactive(data = df[df$type == vName1, ],
      #                            mapping = aes(x = coordinate,y = value,
      #                                          fill = motif,color = motif,
      #                                          tooltip = paste0("Coordinate: ", coordinate, "\nValue: ", round(value, 2)),
      #                                          data_id = coordinate),
      #                            size = 1) +
      #       geom_line_interactive(data = positive_strand_temp,
      #                             mapping = aes(x = coordinate,y = value,
      #                                           linetype = paste0(vName2, " information"),
      #                                           tooltip = paste0("Coordinate: ", coordinate, "\nValue: ", ),
      #                                           data_id = paste0("pos_", coordinate)),
      #                             color = "#868686FF",alpha = alpha) +
      #       geom_line_interactive(data = negative_strand_temp,
      #                             mapping = aes(x = coordinate,y = value,
      #                                           linetype = paste0(vName2, " information"),
      #                                           tooltip = paste0("Coordinate: ", coordinate, "\nValue: ", round(value, 2)),
      #                                           data_id = paste0("neg_", coordinate)),
      #                             color = "#868686FF", alpha = alpha, ) +
      #       labs(linetype = vName2)

      p <- ggplot() +
            geom_col_interactive(data = df[df$type == vName1, ],
                                 mapping = aes(x = coordinate,y = value,
                                               fill = motif,color = motif,
                                               tooltip = paste0("Coordinate: ", coordinate, "\nValue: ", round(value, 2)),
                                               data_id = coordinate),
                                 size = 1) +
            geom_line(data = positive_strand_temp,
                      mapping = aes(x=coordinate,y=value,linetype=paste0(vName2," information")),
                      color = "#868686FF",alpha=alpha) + 
            geom_line(data = negative_strand_temp,
                      mapping = aes(x=coordinate,y=value,linetype=paste0(vName2," information")),
                      color = "#868686FF",alpha=alpha) + 
            labs(linetype = legend_lab_value2)

    }else{
      p <- ggplot(df)+
            geom_col(data = df[df$type==vName1,],mapping = aes(x=coordinate,y=value,fill=motif, color=motif)) +
            labs(fill = legend_lab_motif) +
            geom_line(data = positive_strand_temp,
                      mapping = aes(x=coordinate,y=value,linetype=paste0(vName2," information")),
                      color = "#868686FF",alpha=alpha) + 
            geom_line(data = negative_strand_temp,
                      mapping = aes(x=coordinate,y=value,linetype=paste0(vName2," information")),
                      color = "#868686FF",alpha=alpha) + labs(linetype = legend_lab_value2)
            
    }
    
    # fix order
    p <- p + guides(fill = guide_legend(order = 1),
                    linetype = guide_legend(order = 0),
                    color = "none")


    # reorganize the axis
    p <- p +
      scale_y_continuous(sec.axis = sec_axis(trans = ~rescale(.,c(-value2_max,value2_max)),
                                             name = second_ylab,
                                             breaks = c(-value2_max,
                                                        0,
                                                        value2_max),
                                             labels = c(value2_max,
                                                        0,
                                                        value2_max)),
                         breaks = c((-1)*value1_max, (-0.5)*value1_max, 0, (0.5)*value1_max, value1_max),
                         labels = c(paste0("(3'->5') ", value1_max), 
                                    (0.5)*value1_max, 0, 
                                    (0.5)*value1_max, 
                                    paste0("(5'->3') ",value1_max)))
    # coord_cartesian(ylim = c(-value1_max,value1_max))

  }else{
    # plot the methylation information

    if(interactive){
      p <- ggplot(df) +
            geom_col_interactive(mapping = aes(x = coordinate,y = value,
                                               fill = motif,color = motif,
                                               tooltip = paste0("Coordinate: ", coordinate, "\nValue: ", round(value, 2)),
                                               data_id = coordinate),size = 1)+
            suppressMessages(scale_y_continuous(breaks = c((-1)*value1_max, 
                                                          (-0.5)*value1_max, 
                                                          0, 
                                                          (0.5)*value1_max, 
                                                          value1_max),
                                                labels = c(paste0("(3'->5') ", value1_max), 
                                                          (0.5)*value1_max, 
                                                          0, 
                                                          (0.5)*value1_max, 
                                                          paste0("(5'->3') ",value1_max))))
    }else{
      p <- ggplot(df) +
        geom_col(mapping = aes(x=coordinate,y=value,fill=motif,color=motif))+
        labs(fill = legend_lab_motif) + guides(color = "none") +
        # coord_cartesian(ylim = c(-value1_max,value1_max)) +
        suppressMessages(scale_y_continuous(breaks = c((-1)*value1_max, 
                                                      (-0.5)*value1_max, 
                                                      0, 
                                                      (0.5)*value1_max, 
                                                      value1_max),
                                            labels = c(paste0("(3'->5') ", value1_max), 
                                                      (0.5)*value1_max, 
                                                      0, 
                                                      (0.5)*value1_max, 
                                                      paste0("(5'->3') ",value1_max))))
    }
    
  }

  p <- p + theme_classic()
  

  if(!is.null(highlight)){

    if(n0 ==1){

      if(is(highlight, "list")){

        for(idx in seq_along(highlight)){
            tmp <- highlight[[idx]]
            p <- p + annotate("rect", 
                              xmin = tmp[1], xmax = tmp[2],
                              ymin = -value1_max, ymax = value1_max, 
                              fill = highlight_color, alpha = highlight_alpha)
        }

        }else{
        p <- p + annotate("rect", 
                          xmin = highlight[1], xmax = highlight[2],
                          ymin = -value1_max, ymax = value1_max, 
                          fill = highlight_color, alpha = highlight_alpha)
      }


    }else{

          if(is(highlight, "list")){

        for(idx in seq_along(highlight)){
            tmp <- highlight[[idx]]
            p <- p + annotate("rect", 
                              xmin = tmp[1], xmax = tmp[2],
                              ymin = -value2_max, ymax = value2_max, 
                              fill = highlight_color, alpha = highlight_alpha)
        }

        }else{
        p <- p + annotate("rect", 
                          xmin = highlight[1], xmax = highlight[2],
                          ymin = -value2_max, ymax = value2_max, 
                          fill = highlight_color, alpha = highlight_alpha)
      }

    }

    p <- gginnards::move_layers(x = p,idx = 4L,position = "bottom")

  }

  # facet the plot by sample

  if(length(unique(df$sample)) > 1){
    p <- p + facet_grid(sample~ .) +
          theme(strip.text.y = element_text(angle = angle_of_facet_label,
                                            color = "black",face = "bold",
                                            size = facet_label_text_size))

  }



  if(!strip_border){
    p <- p + theme(strip.background = element_blank())
  }

  # design new x-axis
  p <- p +theme(axis.line.x = element_blank()) +
    geom_hline(yintercept = 0)

  # reorganize the x-axis text
  if(!auto_x_axis){
    p <- p + scale_x_continuous(breaks = c(coordinate[1],
                                           coordinate[floor(length(coordinate)*0.25)],
                                           coordinate[floor(length(coordinate)*0.5)],
                                           coordinate[floor(length(coordinate)*0.75)],
                                           coordinate[length(coordinate)]),
                                labels = c(paste0(df_chr,":",coordinate[1]),
                                           paste0(df_chr,":",coordinate[floor(length(coordinate)*0.25)]),
                                           paste0(df_chr,":",coordinate[floor(length(coordinate)*0.5)]),
                                           paste0(df_chr,":",coordinate[floor(length(coordinate)*0.75)]),
                                           paste0(df_chr,":",coordinate[length(coordinate)])))
  }
  

  # add in label information
  p <- p + labs(title = title, x = xlab, y = ylab) +
    theme(plot.title = element_text(hjust = 0.5,size = title_text_size),
          axis.title = element_text(size = axis_title_text_size,),
          axis.text.x = element_text(vjust = -1.2),
          axis.title.x = element_text(vjust = -1.2),
          axis.title.y.right = element_text(vjust = 5.5))

  # resize text in x and y axis
  p <- p +
    theme(axis.text.y.left = element_text(size = left_y_axis_text_size),
          axis.text.y.right = element_text(size = right_y_axis_text_size),
          axis.text.x = element_text(size = x_axis_text_size))


  # if the color is specified
  if(!is.null(motif_color)){
    # check the length of color
    if(length(motif_color) != length(unique(df$motif))){
      stop("the length of color and the motif should be equal...")
    }
    p <- p + scale_fill_manual(values = motif_color)
  }

  # reorganize the ticks length
  p <- p + theme(axis.ticks.length.y = unit(y_ticks_length, "cm"),
                 axis.ticks.length.x = unit(x_ticks_length, "cm"),
                 plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                 panel.spacing.y = unit(panel_spacing,"cm"))

  # replace the legend
  p <- p + theme(legend.box.spacing = unit(legend_box_spacing,"cm"),
                 legend.position = legend_position)

  if(interactive){
    p <- girafe(ggobj = p, width_svg = width_svg,
                height_svg = height_svg,
                options = list(opts_hover(css = "fill: orange; stroke: black; stroke-width: 2px;"),
                               opts_hover_inv(css = "opacity: 0.5;"),
                               opts_zoom(min = .5, max = 5),
                               opts_tooltip(css = "background-color: white; border: 1px solid black; padding: 5px; border-radius: 3px;")))

  }

  return(p)
}

