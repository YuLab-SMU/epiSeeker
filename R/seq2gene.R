#' Annotate genomic regions to genes in many-to-many mapping
#'
#' This function associates genomic regions with coding genes in a many-to-many mapping. It first maps genomic regions to host genes (either located in exon or intron), proximal genes (located in promoter regions) and flanking genes (located in upstream and downstream within user-specified distance).
#' @title seq2gene
#' @param seq genomic regions in GRanges object
#' @param tssRegion TSS region
#' @param flankDistance flanking search radius
#' @param TxDb TxDb object
#' @param sameStrand logical, whether find nearest/overlap gene in the same strand
#' @return gene vector
#' @export
#' @examples
#' data(seq2gene_result)
#' seq2gene_result
#' @importFrom yulab.utils get_cache_element
#' @importFrom yulab.utils update_cache_item
#' @author Guangchuang Yu
seq2gene <- function(seq, tssRegion, flankDistance, TxDb, sameStrand = FALSE) {
    .epiSeekerEnv(TxDb, item = epiSeekerCache)
    # epiSeekerEnv <- get("epiSeekerEnv", envir=.GlobalEnv)

    ## Exons
    # if ( exists("exonList", envir=epiSeekerEnv, inherits=FALSE) ) {
    #     exonList <- get("exonList", envir=epiSeekerEnv)
    # } else {
    #     exonList <- exonsBy(TxDb)
    #     assign("exonList", exonList, envir=epiSeekerEnv)
    # }
    exonList <- get_cache_element(item = epiSeekerCache, elements = "exonList")
    if (is.null(exonList)) {
        exonList <- exonsBy(TxDb)
        update_cache_item(item = epiSeekerCache, list("exonList" = exonList))
    }
    exons <- getGenomicAnnotation.internal(seq, exonList, type = "Exon", sameStrand = sameStrand)

    ## Introns
    # if ( exists("intronList", envir=epiSeekerEnv, inherits=FALSE) ) {
    #     intronList <- get("intronList", envir=epiSeekerEnv)
    # } else {
    #     intronList <- intronsByTranscript(TxDb)
    #     assign("intronList", intronList, envir=epiSeekerEnv)
    # }
    intronList <- get_cache_element(item = epiSeekerCache, elements = "intronList")

    if (is.null(intronList)) {
        intronList <- intronsByTranscript(TxDb)
        update_cache_item(item = epiSeekerCache, list("intronList" = intronList))
    }

    introns <- getGenomicAnnotation.internal(seq, intronList, type = "Intron", sameStrand = sameStrand)

    genes <- c(exons$gene, introns$gene)
    ## > head(genes)
    ## [1] "uc001aed.3/126789"    "uc001aka.3/440556"    "uc001ako.3/49856"
    ## [4] "uc001alg.3/100133612" "uc009vly.2/390992"    "uc001awv.2/79814"
    genes <- gsub("\\w+\\.*\\d*/(\\d+)", "\\1", genes)
    ## > head(genes)
    ## [1] "126789"    "440556"    "49856"     "100133612" "390992"    "79814"

    features <- getGene(TxDb, by = "gene")
    idx.dist <- getNearestFeatureIndicesAndDistances(seq, features, sameStrand = sameStrand)
    nearestFeatures <- features[idx.dist$index]

    distance <- idx.dist$distance

    nearest_genes <- mcols(nearestFeatures[abs(distance) < flankDistance |
                          (distance > tssRegion[1] & distance < tssRegion[2])])[["gene_id"]]
    genes <- c(genes, nearest_genes)
    return(unique(genes))
}
