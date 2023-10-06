# Create ExperimentHub metadata for cpgea_wgbs_hg38
cpgea_wgbs_hg38_metadata = data.frame(
  Title = "cpgea_wgbs_hg38",
  Description = "A HDF5-backed RangedSummarizedExperiment for WGBS Data 
    (CpG sites only) from matching primary human prostate tumours and 
    normal prostate samples",
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
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "record/8411175/files/se.rds:record/8411175/files/assays.h5"
)
write.csv(x = cpgea_wgbs_hg38_metadata, 
    file = "../extdata/cpgea_wgbs_hg38_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(
    pathToPackage = "~/git_repos/TumourMethData/", 
    fileName = "cpgea_wgbs_hg38_metadata.csv")

# Create ExperimentHub metadata for tcga_wgbs_hg38
tcga_wgbs_hg38_metadata = data.frame(
  Title = "tcga_wgbs_hg38",
  Description = "A HDF5-backed RangedSummarizedExperiment for WGBS Data 
    (CpG sites only) from 39 bladder, breast, colon, glioblastoma, lung, 
    rectal stomach and uterine primary tumour samples and 8 matching normal samples",
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
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "record/8397049/files/se.rds:record/8397049/files/assays.h5"
)
write.csv(tcga_wgbs_hg38_metadata, "../extdata/tcga_wgbs_hg38_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "tcga_wgbs_hg38_metadata.csv")

# Create ExperimentHub metadata for MCRPC
mcrpc_wgbs_hg38_metadata = data.frame(
  Title = "mcrpc_wgbs_hg38",
  Description = "A HDF5-backed RangedSummarizedExperiment for WGBS Data 
    (CpG sites only) from 100 prostate metastases",
  BiocVersion = "3.18",
  Genome = "hg38",
  SourceType = "BED",
  SourceUrl = NA,
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = TRUE,
  DataProvider = "University of California San Francisco",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  Location_Prefix = "https://zenodo.org/",
  RDataClass = "RangedSummarizedExperiment",
  DispatchClass = "H5File",
  RDataPath = "record/8397049/files/se.rds:record/8397049/files/assays.h5"
)
write.csv(mcrpc_wgbs_hg38_metadata, "../extdata/mcrpc_wgbs_hg38_metadata", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "mcrpc_wgbs_hg38_metadata.csv")
