#' TumourMethData: A collection of DNA methylation datasets for human tumour 
#' samples and matching normal samples
#' 
#' TumourMethData collects tumour methylation data from a variety of different 
#' tumour types (and also matching normal samples where available) and produced
#' with different technologies (e.g. WGBS, RRBS and methylation arrays) and 
#' provides them as RangedSummarizedExperiments, facilitating easy extraction 
#' of methylation data for regions of interest. At present, includes the 
#' following datasets: 
#' 
#' * cpgea_wgbs_hg38: WGBS Data from 187 pairs of matching 
#' human prostate tumours and normal prostate samples.
#' * tcga_wgbs_hg38: WGBS Data from 39 bladder, breast, 
#' colon, glioblastoma, lung, rectal stomach and uterine tumour samples and 
#' and 8 matching normal samples.
#' * mcrpc_wgbs_hg38: WGBS data from 100 prostate cancer metastases. 
#' * mcrpc_wgbs_hg38_chr11: Subset of mcrpc_wgbs_hg38 with methylation values
#' for just chromosome 11. 
#' * cao_esophageal_wgbs_hg19: WGBS data for 10 squamous esophageal tumours and 
#' 9 matching normal samples. 
#' @author Richard Heery
#' @docType package
#' @name TumourMethData-package
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import ExperimentHub
#' @import HDF5Array
#' @import rhdf5
#' @md
"_PACKAGE"

#' TumourMethDatasets
#'
#' A table describing the datasets available from TumourMethData.
#'
#'@format A data.frame with one row for each dataset
"TumourMethDatasets"