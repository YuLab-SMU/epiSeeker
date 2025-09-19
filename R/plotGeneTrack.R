#' plot gene track to plot 
#' 
#' @param txdb txdb object, providing gene annotation.
#' @param chr chromosome id.
#' @param start_pos start coordiante of windows.
#' @param end_pos end coordiante of windows.
#' @param xlab x lab.
#' @param ylab y lab.
#' @param x_text_size the size of x text.
#' @param y_text_size the size of y text.
#' @param select_gene show all gene or specifc gene. (1)"all", show all genes. (2) gene symbol, e.g. c("SKAP1", "EFCAB13"). (3) gene id, e.g. c(4831, 55316)
#' @param palette palette, default "Set3".
#' @param fromType from which type of gene name to change gene id. Default: ENTREZID. See ?clusterProfiler::bitr
#' @param highlight a region or a list of region to highlight.
#' @param highlight_color colors of highlight rect. Default "#c6c3c3"
#' @param highlight_alpha alpha of highlight rect.
#' @param OrgDb OrgDb for change gene id to gene symbol.
#' @param show_legend show legend or not.
#' @param auto_x_axis use auto x axis or not.
#' @return ggplot object
#' @importFrom IRanges IRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicFeatures genes
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 scale_fill_brewer
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 coord_cartesian
#' @importFrom rlang sym
#' @importFrom rlang check_installed
#' @examples 
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' plotGeneTrack(txdb = txdb, chr = "chr8", start_pos = 126712193, end_pos = 126713193)
#' @export 
plotGeneTrack <- function(txdb, chr, start_pos, end_pos, xlab = "", ylab = "",
                          x_text_size = 10, y_text_size = 10,
                          select_gene = "all", palette = NULL,  fromType = "ENTREZID",
                          highlight = NULL, highlight_color = "#c6c3c3", highlight_alpha = 0.2,
                          OrgDb = NULL, show_legend = FALSE, auto_x_axis = TRUE) {
    

    # Get genes in the region
    win <- GRanges(seqnames = chr,
                   ranges = IRanges::IRanges(start = start_pos, end = end_pos))
    
    # Get genes related to region
    all_genes <- suppressMessages(GenomicFeatures::genes(txdb))
    gene_df <- data.frame(subsetByOverlaps(x = all_genes, ranges = win, type = "any"))
    gene_df$gene_id <- factor(gene_df$gene_id, levels = unique(gene_df$gene_id))
    
    # Process gene data
    gene_df <- gene_df[,c("seqnames", "gene_id", "start", "end", "strand")]
    colnames(gene_df)[1] <- "chromosome"
    gene_df$forward <- ifelse(gene_df$strand=="+", TRUE, FALSE)

    gene_df$start <- pmax(gene_df$start, start_pos)
    gene_df$end <- pmin(gene_df$end, end_pos)  
    
    # Convert gene IDs
    if(!is.null(OrgDb)){
        # check clusterProfiler install or not 
        rlang::check_installed('clusterProfiler', reason = 'For coverting gene ids.')

        if (requireNamespace("clusterProfiler", quietly = TRUE)){
            changeid <- clusterProfiler::bitr(geneID = gene_df$gene_id, 
                                              fromType = fromType, toType = "SYMBOL",
                                              OrgDb = OrgDb)
        }
        
        colnames(changeid)[1] <- c("gene_id") 
        gene_df <- merge(gene_df, changeid, all.x = TRUE)
        colnames(gene_df)[ncol(gene_df)] <- "gene_symbol"
        legend_label <- "gene_symbol"
        
    }else{

        legend_label <- "gene_id"
    }

    # Subset gene_df
    if(all(select_gene != "all")){

        # gene symbol
        if(!is.null(OrgDb)){
            flag_symbol <- sum(select_gene %in% gene_df$gene_symbol)
        }else{
            flag_symbol <- 0
        }

        # gene id 
        flag_id <- sum(select_gene %in% gene_df$gene_id)

        if(flag_id == length(select_gene)){
            gene_df <- gene_df[gene_df$gene_id %in% select_gene,]
        }else if(flag_symbol == length(select_gene)){
            gene_df <- gene_df[gene_df$gene_symbol %in% select_gene,]
        }else{
            cat("There is no gene selected. Show all genes by default...")
        }

    }

    forward <- NULL

    rlang::check_installed('gggenes', reason = 'For ploting genes structure.')

    if(requireNamespace("gggenes", quietly = TRUE)){

        p <- ggplot(gene_df, aes(xmin = start, xmax = end, y = !!rlang::sym(legend_label), 
                             forward=forward, fill = !!rlang::sym(legend_label))) + 
                gggenes::geom_gene_arrow() +
                coord_cartesian(xlim = c(start_pos, end_pos)) +
                xlim(c(start_pos, end_pos)) +
                gggenes::theme_genes() +
                labs(fill = gsub("_", " ",legend_label), x = xlab, y = ylab) +
                theme(panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    panel.border = element_blank(),
                    plot.title = element_text(hjust = 0.5),
                    axis.text.x = element_text(size = x_text_size),
                    axis.text.y = element_text(size = y_text_size))

    }

    
    
    if(!auto_x_axis){
        p <- p  + scale_x_continuous(breaks = round(as.numeric(quantile(seq(start_pos, end_pos), c(0,0.25,0.5,0.75,1)))),
                                     labels = round(as.numeric(quantile(seq(start_pos, end_pos), c(0,0.25,0.5,0.75,1)))))
    }

    if(!is.null(palette)){
        p <- p + scale_fill_brewer(palette = palette)
    }

    if (!show_legend) {
        p <- p + theme(legend.position = "none")
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

    
    return(p)
}