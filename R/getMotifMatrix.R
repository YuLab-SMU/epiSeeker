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
#' if(FALSE){
#'   require(BSgenome.Hsapiens.UCSC.hg38)
#'   data(pwm_obj)
#'   
#'   region_oi <- GRanges(seqnames = "chr22", 
#'                        ranges = IRanges(start = 10525891, end = 10525991))
#'   motifMatrix <- getMotifMatrix(region = region_oi, 
#'                                 pwm = pwm_obj, 
#'                                 ref_obj = BSgenome.Hsapiens.UCSC.hg38)
#' }
#' 
#' @export 
getMotifMatrix <- function(region, pwm, ref_obj, by = "name"){

    rlang::check_installed('motifmatchr', reason = 'Matching motif...')

    by <- match.arg(by, c("name", "ID"))

    if(by == "name"){
        # Vectorized approach for PWM name processing
        pwm_name <- sapply(pwm, function(x) x@name)
        
        # Handle duplicate PWM names with vectorized approach using ave()
            dup_mask <- duplicated(pwm_name)
            if(any(dup_mask)){
                # Create sequence numbers within each group
                seq_nums <- ave(pwm_name, pwm_name, FUN = seq_along)
                # Only rename duplicates
                pwm_name[dup_mask] <- paste0(pwm_name[dup_mask], "_", seq_nums[dup_mask])
            }
        
        names(pwm) <- pwm_name
    }

    # get start and end position
    start_pos <- start(region)
    end_pos <- end(region)
    chr_name <- as.character(seqnames(region))

    # get seq base
    rlang::check_installed('BSgenome', reason = 'For motif analysis.')

    if (requireNamespace("BSgenome", quietly = TRUE)){
        regionSeqs <- BSgenome::getSeq(ref_obj, region)
    }
    names(regionSeqs) <- as.character(region)

    # get motif position
    if(requireNamespace("motifmatchr", quietly = TRUE)){
        motif_positions <- tryCatch({
            motifmatchr::matchMotifs(pwm, regionSeqs, out = "positions")
        }, error = function(e) {
            return(NULL)
        })
        if(is.null(motif_positions)){
            return(motif_positions)
        }
    }
    
    # Vectorized approach for motif filtering
    motif_lengths <- sapply(motif_positions, function(x) length(x[[1]]))
    sub_list <- motif_positions[motif_lengths > 0]
    
    if(length(sub_list) == 0){
        message("There is no motif match...")
        return(NULL)
    }

    # Pre-allocate list for better performance
    all_list <- vector("list", length(sub_list))
    motif_names <- names(sub_list)
    
    for(i in seq_along(sub_list)){
        tmp_IRange <- sub_list[[i]][[1]]
        mcols(tmp_IRange)[,"motif"] <- motif_names[i]
        
        if(length(tmp_IRange) > 1){
            # Find maximum score more efficiently
            scores <- mcols(tmp_IRange)[,"score"]
            max_score <- max(scores)
            score_idx <- which(scores == max_score)[1]  # Take first occurrence
            all_list[[i]] <- tmp_IRange[score_idx]
        } else {
            all_list[[i]] <- tmp_IRange
        }
    }
    
    # Combine all ranges efficiently
    all_IRange <- as.data.frame(do.call(c, all_list))
    all_IRange$start <- all_IRange$start + start_pos - 1 
    all_IRange$end <- all_IRange$end + start_pos - 1 

    # Vectorized approach for creating motif data frame
    motif_range <- vector("list", nrow(all_IRange))
    
    for(i in seq_len(nrow(all_IRange))){
        start_coord <- all_IRange[i, "start"]
        end_coord <- all_IRange[i, "end"]
        coord_length <- end_coord - start_coord + 1
        
        # Create coordinate sequence efficiently
        tmp_df <- data.frame(
            chr = chr_name,
            coordinate = seq.int(start_coord, end_coord),
            score = all_IRange[i, "score"],
            strand = all_IRange[i, "strand"],
            motif = all_IRange[i, "motif"]
        )
        
        # Adjust score for negative strand
        if(all_IRange[i, "strand"] == "-"){
            tmp_df$score <- (-1) * tmp_df$score
        }
        
        motif_range[[i]] <- tmp_df
    }

    # Efficient rbind using data.table if available
    if(requireNamespace("data.table", quietly = TRUE)){
        motif_df <- data.table::rbindlist(motif_range)
    } else {
        motif_df <- yulab.utils::rbindlist(motif_range)
    }

    attr(motif_df, "range") <- c(start_pos, end_pos)

    return(motif_df)
}