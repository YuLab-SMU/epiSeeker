# Prepare pancancer_atac.rds from the source TCGA pan-cancer ATAC-seq data.
#
# The data originate from:
# Corces MR, Granja JM, Shams S, et al. The chromatin accessibility landscape
# of primary human cancers. Science. 2018;362(6413):eaav1898.
# doi:10.1126/science.aav1898
#
# Source dataset: TCGA-ATAC_PanCan_Log2Norm_Counts.rds, available from:
# https://gdc.cancer.gov/about-data/publications/ATACseq-AWG
# This script generates pancancer_atac.rds from that log2-normalized matrix
# by selecting the target region and averaging signals within each cancer type.
#
# Input: TCGA-ATAC_PanCan_Log2Norm_Counts.rds, a table with seven annotation
# columns, including seqnames, start, and end, followed by sample columns
# containing log2-normalized signals.
#
# Keep only peaks fully contained in chr8:126712193-128412193.
# Derive cancer types from the first four characters of each sample column
# name, removing lowercase x characters from those prefixes, then sort them.
# For each peak, calculate the arithmetic mean of log2-normalized signals
# across samples of each cancer type. Missing values are not removed
# (rowMeans uses na.rm = FALSE); this is not a log2 transform of mean counts.
#
# Output: pancancer_atac.rds, a named list of GRanges objects with strand "*".
# The example dataset contains 23 cancer types with the same 1,555 peak
# intervals in every object. V4 is a cancer-specific peak identifier;
# V5 is the mean signal for that cancer type.

chr <- "chr8"
start_pos <- 126712193
end_pos <- 128412193

atac_data <- readRDS("./raw_data/TCGA-ATAC_PanCan_Log2Norm_Counts.rds")
atac_data$seqnames <- as.character(atac_data$seqnames)
sub_data <- atac_data[
    atac_data$seqnames == chr &
        atac_data$start >= start_pos & atac_data$end <= end_pos,
    , drop = FALSE
]

sample_signals <- sub_data[, 8:ncol(sub_data), drop = FALSE]
cancer_type <- substr(colnames(sample_signals), start = 1, stop = 4)
cancer_type <- gsub("x", "", cancer_type)
cancer_types <- sort(unique(cancer_type))

# Average log2-normalized signals within each cancer type for each peak.
cancer_type_mt <- do.call(cbind, lapply(cancer_types, function(cancer) {
    rowMeans(sample_signals[, cancer_type == cancer, drop = FALSE])
}))
colnames(cancer_type_mt) <- cancer_types

peak_list <- setNames(vector("list", length(cancer_types)), cancer_types)
for (cancer in cancer_types) {
    peak_list[[cancer]] <- GenomicRanges::GRanges(
        seqnames = sub_data$seqnames,
        ranges = IRanges::IRanges(start = sub_data$start, end = sub_data$end),
        strand = "*",
        V4 = paste0(cancer, "_", seq_len(nrow(sub_data))),
        V5 = cancer_type_mt[, cancer]
    )
}

dir.create("./data", recursive = TRUE, showWarnings = FALSE)
saveRDS(peak_list, "./data/pancancer_atac.rds")
