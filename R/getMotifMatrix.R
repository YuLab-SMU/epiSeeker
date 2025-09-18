#' get the information of motif in a range
#' 
#' @param region region object in granges.
#' @param pwm PFMatrixList.
#' @param ref_obj seq reference object. e.g. BSgenome object.
#' @param by show the motif by name or ID.
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom BSgenome getSeq
#' @importFrom GenomeInfoDb seqnames
#' @importFrom rlang check_installed
#' @export 
getMotifMatrix <- function(region, pwm, ref_obj, by = "name"){

    rlang::check_installed('motifmatchr', reason = 'Matching motif...')

    by <- match.arg(by, c("name", "ID"))

    if(by == "name"){
        pwm_name <- vector(mode = "character", length = length(pwm))
        for(i in  seq_len(length(pwm))){
            pwm_name[i] <- pwm[[i]]@name
        }

        new_name <- pwm_name
        for(i in pwm_name[duplicated(pwm_name)]){
            tmp <- pwm_name[pwm_name == i]
            tmp <- paste0(tmp, "_", seq_along(tmp))
            new_name[new_name == i] <- tmp
        }

        names(pwm) <- new_name
    }

    # get start and end position
    start_pos <- start(region)
    end_pos <- end(region)
    chr_name <- as.character(seqnames(region))

    # get seq base
    regionSeqs <- getSeq(ref_obj, region)
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
    
    motif_idx <- vector(length = length(motif_positions))
    for(i in seq_len(length(motif_positions))){
        tmp_idx <- names(motif_positions)[i]
        if(length(motif_positions[[i]][[1]]) == 0){
            motif_idx[i] <- FALSE
        }else{
            motif_idx[i] <- TRUE
        }
    }
    sub_list <- motif_positions[motif_idx]
    if(length(sub_list) == 0){
        cat("There is no motif match...")
        return(NULL)
    }

    for(i in seq_len(length(sub_list))){
        mcols(sub_list[[i]][[1]])[,"motif"] <- names(sub_list)[i]
    }

    all_list <- list()
    for(i in seq_len(length(sub_list))){
        tmp_IRange <- sub_list[[i]][[1]]
        if(length(tmp_IRange) >1){
            score_idx <- which(mcols(tmp_IRange)[,"score"] == max(mcols(tmp_IRange)[,'score']))
            if(length(score_idx) > 1){
                score_idx <- score_idx[1]
            }
            all_list[[i]] <- tmp_IRange[score_idx]
        }else{
            all_list[[i]] <- tmp_IRange
        }
        
    }
    all_IRange <- as.data.frame(do.call(c,all_list))
    all_IRange$start <- all_IRange$start + start_pos - 1 
    all_IRange$end <- all_IRange$end + start_pos - 1 

    # create a data frame
    motif_range <- list()

    for(i in seq_len(nrow(all_IRange))){
        # tmp_bak_df <- data.frame(chr = chr_name,
        #                          coordinate = seq(from = start_pos,
        #                                           to = end_pos,
        #                                           by = 1),
        #                          score = 0,
        #                          strand = all_IRange[i, "strand"],
        #                          motif = all_IRange[i,"motif"])
            

        tmp_df <- data.frame(chr = chr_name,
                             coordinate = seq(from = all_IRange[i, "start"],
                                              to = all_IRange[i, "end"],
                                              by = 1))
        tmp_df$score <- all_IRange[i, "score"]
        tmp_df$strand <- all_IRange[i, "strand"]
        tmp_df$motif <- all_IRange[i,"motif"]
        if(all_IRange[i, "strand"] == "-"){
            tmp_df$score <- (-1)*tmp_df$score
        }
        
        # tmp_match_indices <- match(tmp_df$coordinate, tmp_bak_df$coordinate)
        # tmp_valid_indices <- which(!is.na(tmp_match_indices))
        # tmp_bak_df[tmp_match_indices[tmp_valid_indices], ] <- tmp_df[tmp_valid_indices, ]
        # motif_range[[i]] <- tmp_bak_df

        motif_range[[i]] <- tmp_df
    }

    motif_df <- do.call(rbind, motif_range)

    attr(motif_df, "range") <- c(start_pos, end_pos)

    return(motif_df)
}