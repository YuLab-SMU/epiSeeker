# epiSeeker 0.99.10

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
