#' get the information of motif in a range
#' 
#' @param region region object in granges.
#' @param pwm PFMatrixList.
#' @param ref_obj seq reference object. e.g. BSgenome object.
#' @param by show the motif by name or ID.
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom rlang check_installed
#' @return score matrix
#' @examples 
#' require(BSgenome.Hsapiens.UCSC.hg38)
#' data(pwm_obj)
#'   
#' region_oi <- GRanges(seqnames = "chr22", 
#'                      ranges = IRanges(start = 10525891, end = 10525991))
#' motifMatrix <- getMotifMatrix(region = region_oi, 
#'                               pwm = pwm_obj[c(45,120,170)], 
#'                               ref_obj = BSgenome.Hsapiens.UCSC.hg38)
#' 
#' @export 
getMotifMatrix <- function(region, 
                           pwm, 
                           ref_obj, 
                           by = "name"){

    rlang::check_installed(c('BSgenome', 'universalmotif'), reason = 'Matching motif...')

    by <- match.arg(by, c("name", "ID"))
    
    # get start and end position
    start_pos <- start(region)
    end_pos <- end(region)
    chr_name <- as.character(seqnames(region))

    # get seq base
    if (requireNamespace("BSgenome", quietly = TRUE)){
        regionSeqs <- BSgenome::getSeq(ref_obj, region)
    }

    # get motif position
    if(requireNamespace("universalmotif", quietly = TRUE)){
        motif_df <- as.data.frame(universalmotif::scan_sequences(pwm, regionSeqs, respect.strand = TRUE))
    }
    
    motif_df$len <- motif_df$stop - motif_df$start
    motif_df$start <- motif_df$start + start_pos - 1 
    motif_df$stop <- motif_df$stop + start_pos - 1 

    if(by == "name"){
        motif_df$index <- motif_df$motif
    }else{
        motif_df$index <- names(pwm)[motif_df$motif.i]
    }
    
    motif_list <- lapply(seq_len(nrow(motif_df)), function(i){

        s  <- motif_df$start[i]
        e  <- motif_df$stop[i]
        sc <- motif_df$score[i]
        st <- motif_df$strand[i]
        mt <- motif_df$index[i]
        seg_id <- i

        coords <- seq.int(s, e)

        df <- data.frame(
            chr        = chr_name,
            coordinate = coords,
            score      = sc,
            strand     = st,
            motif      = mt,
            seg_id     = seg_id
        )
        if(st == "-"){
            df$score <- -df$score
        }

        df
    })

    # Efficient rbind using data.table if available
    if(requireNamespace("data.table", quietly = TRUE)){
        result <- data.table::rbindlist(motif_list)
    } else {
        result <- yulab.utils::rbindlist(motif_list)
    }

    attr(result, "range") <- c(start_pos, end_pos)

    return(result)
}