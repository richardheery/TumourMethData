TumourMethDatasets = data.frame(
  dataset_name = c("cpgea_wgbs_hg38", "tcga_wgbs_hg38"), 
  technology = c("WGBS", "WGBS"),
  genome_build = c("hg38", "hg38"),
  tumour_site = c("prostate", "various"),
  number_tumour_samples = c(187, 39), 
  number_normal_samples = c(187, 8),
  coverage_available = c(FALSE, FALSE),
  notes = c("", ""), 
  original_publication = c(
    "A genomic and epigenomic atlas of prostate cancer in Asian populations; Nature; 2020", 
    "DNA methylation loss in late-replicating domains is linked to mitotic cell division; Nature genetics; 2018"
  )
)
usethis::use_data(TumourMethDatasets, overwrite = T)