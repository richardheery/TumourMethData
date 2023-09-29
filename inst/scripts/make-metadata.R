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
write.csv(cpgea_wgbs_hg38_metadata, "../extdata/metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata("~/git_repos/TumourMethData/")
