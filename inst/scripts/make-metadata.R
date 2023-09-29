# Create ExperimentHub metadata for cpgea_wgbs_hg38
cpgea_wgbs_hg38_metadata = data.frame(
  Title = "cpgea_wgbs_hg38",
  Description = "A HDF5-backed RangedSummarizedExperiment for WGBS Data (CpG sites only) from matching human prostate tumours and normal prostate samples",
  BiocVersion = "3.18",
  Genome = "hg38",
  SourceType = "BED",
  SourceUrl = "https://wangftp.wustl.edu/~hlee/SMMU/PC/WGBS_bedGraph/",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = TRUE,
  DataProvider = "Washington University School of Medicine",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "RangedSummarizedExperiment",
  DispatchClass = "H5File",
  RDataPath = "TumourMethData"
)
write.csv(cpgea_wgbs_hg38_metadata, "../extdata/cpgea_wgbs_hg38_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "cpgea_wgbs_hg38_metadata.csv")

# Create ExperimentHub metadata for tcga_wgbs_hg38
tcga_wgbs_hg38_metadata = data.frame(
  Title = "tcga_wgbs_hg38",
  Description = "A HDF5-backed RangedSummarizedExperiment for WGBS Data (CpG sites only) from 39 bladder, breast, colon, glioblastoma, lung, rectal stomach and uterine tumour samples and 8 matching normal samples",
  BiocVersion = "3.18",
  Genome = "hg38",
  SourceType = "BED",
  SourceUrl = "https://zwdzwd.s3.amazonaws.com/directory_listing/trackHubs_TCGA_WGBS_hg38.html?prefix=trackHubs/TCGA_WGBS/hg38/bed/",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = TRUE,
  DataProvider = "Center for Epigenetics, Van Andel Research Institute, Grand Rapids, MI, USA",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "RangedSummarizedExperiment",
  DispatchClass = "H5File",
  RDataPath = "TumourMethData"
)
write.csv(tcga_wgbs_hg38_metadata, "../extdata/tcga_wgbs_hg38_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "tcga_wgbs_hg38_metadata.csv")
