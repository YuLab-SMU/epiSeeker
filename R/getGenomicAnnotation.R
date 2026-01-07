updateGenomicAnnotation <- function(peaks, genomicRegion, type, anno, sameStrand = FALSE) {
    hits <- getGenomicAnnotation.internal(peaks, genomicRegion, type, sameStrand = sameStrand)
    if (length(hits) > 1) {
        hitIndex <- hits$queryIndex
        anno[["annotation"]][hitIndex] <- hits$annotation
        anno[["detailGenomicAnnotation"]][hitIndex, type] <- TRUE
    }
    return(anno)
}


#' Get Genomic Annotation of peaks
#'
#' @title getGenomicAnnotation
#' @param peaks peaks in GRanges object
#' @param distance distance of peak to TSS
#' @param tssRegion tssRegion, default is -3kb to +3kb
#' @param TxDb TxDb object
#' @param level one of gene or transcript
#' @param genomicAnnotationPriority genomic Annotation Priority
#' @param sameStrand whether annotate gene in same strand
#' @importFrom GenomicFeatures threeUTRsByTranscript
#' @importFrom GenomicFeatures fiveUTRsByTranscript
#' @importFrom BiocGenerics unstrand
#' @importFrom yulab.utils get_cache_element
#' @importFrom yulab.utils update_cache_item
#' @return character vector
#' @author G Yu
getGenomicAnnotation <- function(peaks,
                                 distance,
                                 tssRegion = c(-3000, 3000),
                                 TxDb,
                                 level,
                                 genomicAnnotationPriority,
                                 sameStrand = FALSE) {
    .epiSeekerEnv(TxDb, item = epiSeekerCache)

    anno <- .init_anno(length(distance))

    genomicAnnotationPriority <- rev(genomicAnnotationPriority)
    for (AP in genomicAnnotationPriority) {
        anno <- .update_anno_by_priority(AP, peaks, distance, tssRegion, TxDb, anno, sameStrand)
    }

    anno <- .finalize_genic_anno(anno)
    annotation <- anno[["annotation"]]
    detailGenomicAnnotation <- anno[["detailGenomicAnnotation"]]

    features <- getGene(TxDb, by = level)
    dd2 <- .get_distance_to_gene_end(peaks, features, sameStrand)

    dsd <- getOption("epiSeeker.downstreamDistance") %||% 3000

    annotation <- .add_downstream_anno(annotation, dd2, dsd)
    annotation[which(annotation == "Intergenic")] <- "Distal Intergenic"

    downstreamIndex <- dd2 > 0 & dd2 < dsd
    detailGenomicAnnotation[downstreamIndex, "downstream"] <- TRUE
    detailGenomicAnnotation[which(annotation == "Distal Intergenic"), "distal_intergenic"] <- TRUE

    return(list(annotation = annotation, detailGenomicAnnotation = detailGenomicAnnotation))
}

.init_anno <- function(len) {
    annotation <- rep(NA, len)
    flag <- rep(FALSE, len)
    detailGenomicAnnotation <- data.frame(
        genic = flag,
        Intergenic = flag,
        Promoter = flag,
        fiveUTR = flag,
        threeUTR = flag,
        Exon = flag,
        Intron = flag,
        downstream = flag,
        distal_intergenic = flag
    )
    list(annotation = annotation, detailGenomicAnnotation = detailGenomicAnnotation)
}

.update_anno_by_priority <- function(AP, peaks, distance, tssRegion, TxDb, anno, sameStrand) {
    if (AP == "Intron") {
        intronList <- get_intronList(item = epiSeekerCache)
        anno <- updateGenomicAnnotation(peaks, intronList, "Intron", anno, sameStrand = sameStrand)
    } else if (AP == "Exon") {
        exonList <- get_exonList(item = epiSeekerCache)
        anno <- updateGenomicAnnotation(peaks, exonList, "Exon", anno, sameStrand = sameStrand)
    } else if (AP == "3UTR") {
        threeUTRList <- get_cache_element(item = epiSeekerCache, elements = "threeUTRList")
        if (is.null(threeUTRList)) {
            threeUTRList <- threeUTRsByTranscript(TxDb)
            update_cache_item(item = epiSeekerCache, list("threeUTRList" = threeUTRList))
        }
        anno <- updateGenomicAnnotation(peaks, threeUTRList, "threeUTR", anno, sameStrand = sameStrand)
    } else if (AP == "5UTR") {
        fiveUTRList <- get_cache_element(item = epiSeekerCache, elements = "fiveUTRList")
        if (is.null(fiveUTRList)) {
            fiveUTRList <- fiveUTRsByTranscript(TxDb)
            update_cache_item(item = epiSeekerCache, list("fiveUTRList" = fiveUTRList))
        }
        anno <- updateGenomicAnnotation(peaks, fiveUTRList, "fiveUTR", anno, sameStrand = sameStrand)
    } else if (AP == "Promoter") {
        anno <- .update_promoter_anno(anno, distance, tssRegion)
    } else {
        anno[["annotation"]][is.na(anno[["annotation"]])] <- "Intergenic"
    }
    return(anno)
}

.update_promoter_anno <- function(anno, distance, tssRegion) {
    annotation <- anno[["annotation"]]
    tssIndex <- distance >= tssRegion[1] & distance <= tssRegion[2]
    annotation[tssIndex] <- "Promoter"
    anno$detailGenomicAnnotation[tssIndex, "Promoter"] <- TRUE

    pm <- max(abs(tssRegion))
    if (pm / 1000 >= 2) {
        dd <- seq_len(ceiling(pm / 1000)) * 1000
        for (i in seq_len(length(dd))) {
            if (i == 1) {
                lbs <- paste("Promoter", " (<=", dd[1] / 1000, "kb)", sep = "")
                annotation[abs(distance) <= dd[1] & annotation == "Promoter"] <- lbs
            } else {
                lbs <- paste("Promoter", " (", dd[i - 1] / 1000, "-", dd[i] / 1000, "kb)", sep = "")
                annotation[abs(distance) <= dd[i] & abs(distance) > dd[i - 1] & annotation == "Promoter"] <- lbs
            }
        }
    }
    anno[["annotation"]] <- annotation
    return(anno)
}

.finalize_genic_anno <- function(anno) {
    detailGenomicAnnotation <- anno[["detailGenomicAnnotation"]]
    genicIndex <- which(apply(detailGenomicAnnotation[, c("Exon", "Intron")], 1, any))
    detailGenomicAnnotation[-genicIndex, "Intergenic"] <- TRUE
    detailGenomicAnnotation[genicIndex, "genic"] <- TRUE
    anno[["detailGenomicAnnotation"]] <- detailGenomicAnnotation
    return(anno)
}

.get_distance_to_gene_end <- function(peaks, features, sameStrand) {
    if (sameStrand) {
        idx <- follow(peaks, features)
    } else {
        idx <- follow(peaks, BiocGenerics::unstrand(features))
    }

    na.idx <- which(is.na(idx))
    if (length(na.idx)) {
        idx_no_na <- idx[-na.idx]
        peaks_no_na <- peaks[-na.idx]
        peF <- features[idx_no_na]
        dd <- ifelse(strand(peF) == "+", start(peaks_no_na) - end(peF), end(peaks_no_na) - start(peF))
        dd2 <- numeric(length(idx))
        dd2[-na.idx] <- dd
    } else {
        peF <- features[idx]
        dd2 <- ifelse(strand(peF) == "+", start(peaks) - end(peF), end(peaks) - start(peF))
    }
    return(dd2)
}

.add_downstream_anno <- function(annotation, dd2, dsd) {
    if (dsd / 1000 <= 1) {
        j <- which(annotation == "Intergenic" & abs(dd2) <= dsd & dd2 != 0)
        if (length(j) > 0) {
            annotation[j] <- paste("Downstream (<=", dsd, "bp)", sep = "")
        }
    } else {
        for (i in seq_len((dsd / 1000))) {
            j <- which(annotation == "Intergenic" & abs(dd2) <= i * 1000 & dd2 != 0)
            if (length(j) > 0) {
                lbs <- if (i == 1) "Downstream (<1kb)" else paste("Downstream (", i - 1, "-", i, "kb)", sep = "")
                annotation[j] <- lbs
            }
        }
        z <- which(annotation == "Intergenic" & abs(dd2) <= dsd & dd2 != 0)
        if (length(z) > 0) {
            annotation[z] <- paste("Downstream (", dsd / 1000, "kb-", dsd, "bp)", sep = "")
        }
    }
    return(annotation)
}


#' @import IRanges
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @importFrom BiocGenerics unstrand
getGenomicAnnotation.internal <- function(peaks, genomicRegion, type, sameStrand = FALSE) {
    GRegion <- unlist(genomicRegion)
    GRegionLen <- elementNROWS(genomicRegion)

    names(GRegionLen) <- names(genomicRegion)
    GRegion$gene_id <- rep(names(genomicRegion), times = GRegionLen)


    if (type == "Intron") {
        gr2 <- GRegion[!duplicated(GRegion$gene_id)]
        strd <- as.character(strand(gr2))
        len <- GRegionLen[GRegionLen != 0]

        GRegion$intron_rank <- lapply(seq_along(strd), function(i) {
            rank <- seq(1, len[i])
            if (strd[i] == "-") {
                rank <- rev(rank)
            }
            return(rank)
        }) %>% unlist()
    }

    if (type == "Intron" || type == "Exon") {
        nn <- TXID2EG(names(genomicRegion))
        names(GRegionLen) <- nn
        GRegion$gene_id <- rep(nn, times = GRegionLen)
    }

    ## find overlap
    if (sameStrand) {
        GRegionHit <- findOverlaps(peaks, GRegion)
    } else {
        GRegionHit <- findOverlaps(peaks, BiocGenerics::unstrand(GRegion))
    }

    if (length(GRegionHit) == 0) {
        return(NA)
    }
    qh <- queryHits(GRegionHit)
    hit.idx <- getFirstHitIndex(qh)
    GRegionHit <- GRegionHit[hit.idx]
    queryIndex <- queryHits(GRegionHit)
    subjectIndex <- subjectHits(GRegionHit)

    hits <- GRegion[subjectIndex]
    geneID <- hits$gene_id

    if (type == "Intron") {
        anno <- paste(type, " (", geneID, ", intron ", hits$intron_rank,
            " of ", GRegionLen[geneID], ")",
            sep = ""
        )
    } else if (type == "Exon") {
        anno <- paste(type, " (", geneID, ", exon ", hits$exon_rank,
            " of ", GRegionLen[geneID], ")",
            sep = ""
        )
    } else if (type == "fiveUTR") {
        anno <- "5' UTR"
    } else if (type == "threeUTR") {
        anno <- "3' UTR"
    } else {
        anno <- type
    }
    res <- list(queryIndex = queryIndex, annotation = anno, gene = geneID)
    return(res)
}
