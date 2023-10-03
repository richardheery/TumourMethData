#' TumourMethData: A collection of DNA methylation datasets for human tumour 
#' samples and matching normal samples
#' 
#' TumourMethData collects tumour methylation data from a variety of different 
#' tumour types (and also matching normal samples where available) and produced
#' with different technologies (e.g. WGBS, RRBS and methylation arrays) and 
#' provides them as RangedSummarizedExperiments, facilitating easy extraction 
#' ofmethylation data for regions of interest. At present, includes the 
#' following datasets: 
#' 
#' * cpgea_wgbs_hg38 : WGBS Data (CpG sites only) from 187 pairs of matching 
#' human prostate tumours and normal prostate samples from Li, Jing, et al. 
#' "A genomic and epigenomic atlas of prostate cancer in Asian populations." 
#' Nature 580.7801 (2020): 93-99.'
#' * tcga_wgbs_hg38 : WGBS Data (CpG sites only) from 39 bladder, breast, 
#' colon, glioblastoma, lung, rectal stomach and uterine tumour samples and 
#' and 8 matching normal samples from "Zhou, Wanding, et al. "DNA methylation 
#' loss in late-replicating domains is linked to mitotic cell division."
#' Nature genetics c50.4 (2018): 591-602." 
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