# epiSeeker 1.1.3

+ `annotateSeq()` now reports `geneChr` and `geneStrand` as characters instead of
  factor codes. `as.data.frame()` returns 'seqnames'/'strand' as factors and
  assigning a factor into `mcols()` dropped the class, so both columns came out
  as integers (e.g. `geneStrand` = 1/2 instead of +/-, and a wrong `geneChr`
  whenever the seqlevels were not in numeric order). (2026-09-15, Tue)
+ `annotateSeq()` now warns when peaks are dropped because no feature of `TxDb`
  can be found for them (e.g. peaks on contigs/scaffolds without genes, or a
  seqlevels style mismatch). Previously they disappeared silently. (2026-09-15, Tue)
+ `annotateSeq(..., sameStrand = TRUE)` no longer assigns a peak to a feature on
  the opposite strand. Overlap detection in `getNearestFeatureIndicesAndDistances()`
  was calling `findOverlaps()` with `unstrand(features)`, which bypassed `sameStrand`
  and overrode the strand-aware nearest-feature result. Peaks with ambiguous
  strand (`*`) are unaffected and still match features on any strand.
  (2026-09-15, Tue)
+ `plotAnnoBar()` no longer uses the deprecated `ggplot2::aes_string()`. It follows
  the tidy evaluation idiom already used by `plotDistToTSS()`. (2026-09-15, Tue)
+ dropped unused `@importFrom` directives with no call site (`aplot::xlim2`,
  `ggplot2::geom_segment`, `ggplot2::geom_text`, `ggplot2::scale_fill_hue`,
  `utils::getFromNamespace`) and a duplicate `ggplot2::scale_fill_brewer`. The
  `grid::unit()` import, previously declared in `upsetplot()` which never used it,
  now sits in `plotBmProf()` where `unit()` is actually called, keeping `grid` a
  genuine import. (2026-09-17, Thu)

# epiSeeker 1.1.2

+ documentation build migrated to roxygen2 8.1.0 (re-generated Rd pages and
  NAMESPACE; fixed `csAnno` class `@aliases` to be a single line); the exported
  API is unchanged
+ fixed bug in `getNearestFeatureIndicesAndDistances()` where results for
  `overlap=="all"` were silently overridden by the `overlap=="TSS"` branch,
  causing the two modes to behave identically (2026-09-07, Mon)

# epiSeeker 1.1.1

+ make csAnno subset robust to GRanges metadata columns (2027-07-07, Tue)

# epiSeeker 1.0.0

+ Bioconductor RELEASE_3_23 (2026-04-29, Wed)

# epiSeeker 0.99.12

+ add more runnable examples and increase the coverage of tests (44.18%) (2025-11-30, Sun)
+ use new demo data for base modification (`demo_bmdata`) and move some dependency packages from 'Imports' to 'Suggests' (2025-11-10, Mon)
+ consistently use 'hg38' for demo with new demo file, a small subset derived from GSM6418464 (2025-11-07, Fri)
+ fixed Notes reported by `BiocCheck` and add 'GeneRegulation' to 'biocViews' (2025-11-07, Fri)
+ update test for the changes of TxDb.Hsapiens.UCSC.hg19.knownGene (2025-11-06, Thu)
+ new cache mechanism from 'yulab.utils' (2025-10-15, Wed)
+ fixed R check (2025-10-05, Sun)
+ remove dependent of `Vennerable` package (2025-09-27, Sat)
+ epiSeeker inherits from ChIPseeker to support analysis of multi-omics epigenomic data (2025-09-19, Fri)
    - epiSeeker inherits functions from ChIPseeker to plot peak profile, plot coverage, peak annotaions, visualization of annotation, data mining and statistic for peak overla
    - epiSeeker merges multiple functions to plot profile (getTagMatrix + plotPeakProf/plotPeakHeatmap), and add featrues of dealing missing data and different statistic methods to combind data
    - epiSeeker provides new features for plotCov to plot cluster tree and peak co-accessibility
    - epiSeeekr provides interactive function for user to explore data
    - epiSeeker provides functions to plot gene structure, base modification and motif
    - see more on <https://github.com/YuLab-SMU/ChIPseeker/blob/devel/NEWS.md>
