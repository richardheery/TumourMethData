TumourMethDatasets = data.frame(
  dataset_name = c("cpgea_wgbs_hg38", "tcga_wgbs_hg38", "mcrpc_wgbs_hg38", "mcrpc_wgbs_hg38_chr11", "cao_esophageal_wgbs_hg19", "target_rhabdoid_wgbs_hg19"), 
  cancer_type = c("prostate", "various", "prostate", "prostate", "esophageal", "rhabdoid"),
  technology = c("WGBS", "WGBS", "WGBS", "WGBS", "WGBS", "WGBS"),
  genome_build = c("hg38", "hg38", "hg38", "hg38", "hg19", "hg19"),
  number_tumour_samples = c(187, 39, 100, 100, 10, 69), 
  number_normal_samples = c(187, 8, 0, 0, 9, 0),
  wgbs_coverage_available = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
  dataset_size_gb = c(40, 5.4, 16, 0.76, 2, 4.5),
  transcript_counts_available = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), 
  notes = c("", "", "", "This dataset is a subset of the data in mcrpc_wgbs_hg38 for example purposes", "",
    "Methylation values are not as precise as in other datasets. The original 
    methylation values were integers between 0 and 10 with separate values for the C and G positions of each CpG site.
    The mean of the C and G values for each CpG sie was divided by 10 to produce the methylation values here, 
    with CpG sites missing methylation values for either to C or G given an NA value"), 
  original_publication = c(
    "A genomic and epigenomic atlas of prostate cancer in Asian populations; Nature; 2020", 
    "DNA methylation loss in late-replicating domains is linked to mitotic cell division; Nature genetics; 2018",
    "The DNA methylation landscape of advanced prostate cancer; Nature genetics; 2020",
    "The DNA methylation landscape of advanced prostate cancer; Nature genetics; 2020",
    "Multi-faceted epigenetic dysregulation of gene expression promotes esophageal squamous cell carcinoma; Nature communications; 2020",
    "Genome-Wide Profiles of Extra-cranial Malignant Rhabdoid Tumors Reveal Heterogeneity and Dysregulated Developmental Pathways; Cancer Cell; 2016"
  )
)
usethis::use_data(TumourMethDatasets, overwrite = T)

.experimenthub_ids = data.frame(
  wgbs = c("EH8524", "EH8425", "EH8526", "EH8527", "EH8528"),
  rnaseq = c("", "", "", "", ""), 
  row.names = c("tcga_wgbs_hg38", "mcrpc_wgbs_hg38", "cpgea_wgbs_hg38", "cao_esophageal_wgbs_hg19", "mcrpc_wgbs_hg38_chr11") 
)
usethis::use_data(.experimenthub_ids, overwrite = T, internal = T)
