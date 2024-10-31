# Create ExperimentHub metadata for roadmap_wgbs_hg38
roadmap_wgbs_hg38_metadata = data.frame(
  Title = "roadmap_wgbs_hg38",
  Description = "A HDF5-backed RangedSummarizedExperiment for WGBS Data 
    (CpG sites only) for 38 normal human tissue samples from the Roadmap Epigenomics Project",
  BiocVersion = "3.20",
  Genome = "hg38",
  SourceType = "BED",
  SourceUrl = "https://www.encodeproject.org/",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = TRUE,
  DataProvider = "NIH Roadmap Epigenomics Mapping Consortium",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "RangedSummarizedExperiment",
  DispatchClass = "H5File",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "record/13902805/files/se.rds:record/13902805/files/assays.h5"
)
write.csv(roadmap_wgbs_hg38_metadata, "../extdata/roadmap_wgbs_hg38_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "roadmap_wgbs_hg38_metadata.csv")

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

# Create ExperimentHub metadata for MCRPC
mcrpc_wgbs_hg38_metadata = data.frame(
  Title = "mcrpc_wgbs_hg38",
  Description = "A HDF5-backed RangedSummarizedExperiment for WGBS Data 
    (CpG sites only) from 100 castration-resistant prostate cancer metastases",
  BiocVersion = "3.18",
  Genome = "hg38",
  SourceType = "BED",
  SourceUrl = "https://www.box.com/",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = TRUE,
  DataProvider = "University of California San Francisco",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "RangedSummarizedExperiment",
  DispatchClass = "H5File",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "record/8413802/files/se.rds:record/8413802/files/assays.h5"
)
write.csv(mcrpc_wgbs_hg38_metadata, "../extdata/mcrpc_wgbs_hg38_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "mcrpc_wgbs_hg38_metadata.csv")

# Create ExperimentHub metadata for MCRPC
mcrpc_wgbs_hg38_chr11_metadata = data.frame(
  Title = "mcrpc_wgbs_hg38_chr11",
  Description = "A HDF5-backed RangedSummarizedExperiment for WGBS Data for 
  chromosome 11 (CpG sites only) from 100 castration-resistant prostate cancer 
  metastases",
  BiocVersion = "3.18",
  Genome = "hg38",
  SourceType = "BED",
  SourceUrl = "https://www.box.com/",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = TRUE,
  DataProvider = "University of California San Francisco",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "RangedSummarizedExperiment",
  DispatchClass = "H5File",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "record/8432665/files/se.rds:record/8432665/files/assays.h5"
)
write.csv(mcrpc_wgbs_hg38_chr11_metadata, "../extdata/mcrpc_wgbs_hg38_chr11_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "mcrpc_wgbs_hg38_chr11_metadata.csv")

# Create ExperimentHub metadata for Cao Esophageal WGBS data
cao_esophageal_wgbs_hg19_metadata = data.frame(
  Title = "cao_esophageal_wgbs_hg19",
  Description = "A HDF5-backed RangedSummarizedExperiment for WGBS Data 
    (CpG sites only) from 10 esophageal squamous carcinomas and 9 matching 
    normal esophageal samples",
  BiocVersion = "3.18",
  Genome = "hg19",
  SourceType = "BED",
  SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149608",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = TRUE,
  DataProvider = "University of California San Francisco",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "RangedSummarizedExperiment",
  DispatchClass = "H5File",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "record/8414527/files/se.rds:record/8414527/files/assays.h5"
)
write.csv(cao_esophageal_wgbs_hg19_metadata, "../extdata/cao_esophageal_wgbs_hg19_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "cao_esophageal_wgbs_hg19_metadata.csv")

# Create ExperimentHub metadata for TARGET rhabdoid data 
target_rhabdoid_wgbs_hg19_metadata = data.frame(
  Title = "target_rhabdoid_wgbs_hg19",
  Description = "A HDF5-backed RangedSummarizedExperiment for WGBS Data 
    (CpG sites only) from 69 childhood rhabdoid tumours",
  BiocVersion = "3.18",
  Genome = "hg19",
  SourceType = "BigWig",
  SourceUrl = "https://www.globus.org/",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = TRUE,
  DataProvider = "Canada's Michael Smith Genome Sciences Centre",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "RangedSummarizedExperiment",
  DispatchClass = "H5File",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "record/10019178/files/se.rds:record/10019178/files/assays.h5"
)
write.csv(target_rhabdoid_wgbs_hg19_metadata, "../extdata/target_rhabdoid_wgbs_hg19_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "target_rhabdoid_wgbs_hg19_metadata.csv")

# Create ExperimentHub metadata for TCGA array data
tcga_450k_array_hg19_metadata = data.frame(
  Title = "tcga_450k_array_hg19",
  Description = "A HDF5-backed RangedSummarizedExperiment for Infinium HumanMethylation450 BeadChip 
    array methylation data for all TCGA samples",
  BiocVersion = "3.20",
  Genome = "hg19",
  SourceType = "TSV",
  SourceUrl = "https://gdc.cancer.gov/",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = TRUE,
  DataProvider = "Genomic Data Commons",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "RangedSummarizedExperiment",
  DispatchClass = "H5File",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "record/10019178/files/se.rds:record/10019178/files/assays.h5"
)
write.csv(tcga_450k_array_hg19_metadata, "../extdata/tcga_450k_array_hg19_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "tcga_450k_array_hg19_metadata.csv")

### Create metadata files for RNA-seq count datasets

# Create ExperimentHub metadata for tcga_transcript_counts
tcga_wgbs_sample_transcript_counts_metadata = data.frame(
  Title = "tcga_wgbs_sample_transcript_counts",
  Description = "Transcript counts for 23 tumour and 7 normal samples in 
  tcga_wgbs_hg38 quantified using Kallisto with Gencode version 38 
    transcript annotation",
  BiocVersion = "3.18",
  Genome = "hg38",
  SourceType = "FASTQ",
  SourceUrl = "https://portal.gdc.cancer.gov/repository",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = NA,
  DataProvider = "Center for Epigenetics, Van Andel Research Institute, Grand Rapids, MI, USA",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "data.frame",
  DispatchClass = "FilePath",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "record/8422703/files/tcga_wgbs_sample_transcript_counts.tsv.gz"
)
write.csv(tcga_transcript_counts_metadata, "../extdata/tcga_transcript_counts_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "tcga_transcript_counts_metadata.csv")

# Create ExperimentHub metadata for cpgea_transcript_counts
cpgea_transcript_counts_metadata = data.frame(
  Title = "cpgea_transcript_counts",
  Description = "Transcript counts for 126 pairs of the prostate tumour and 
  matching normal prostate samples from cpgea_wgbs_hg38 quantified using 
  Kallisto with Gencode version 38 transcript annotation",
  BiocVersion = "3.18",
  Genome = "hg38",
  SourceType = "FASTQ",
  SourceUrl = "ftp://human.big.ac.cn/",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = NA,
  DataProvider = "Washington University School of Medicine",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "data.frame",
  DispatchClass = "FilePath",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "record/8422703/files/cpgea_transcript_counts.tsv.gz"
)
write.csv(x = cpgea_transcript_counts_metadata, 
    file = "../extdata/cpgea_transcript_counts_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(
    pathToPackage = "~/git_repos/TumourMethData/", 
    fileName = "cpgea_transcript_counts_metadata.csv")

# Create ExperimentHub metadata for MCRPC transcript counts
mcrpc_transcript_counts_metadata = data.frame(
  Title = "mcrpc_transcript_counts",
  Description = "Transcript counts for 99 samples in mcrpc_wgbs_hg38 
  quantified using Kallisto with Gencode version 38 transcript annotation",
  BiocVersion = "3.18",
  Genome = "hg38",
  SourceType = "FASTQ",
  SourceUrl = "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001648.v2.p1",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = NA,
  DataProvider = "University of California San Francisco",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "data.frame",
  DispatchClass = "FilePath",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "record/8422703/files/mcrpc_transcript_counts.tsv.gz"
)
write.csv(mcrpc_transcript_counts_metadata, "../extdata/mcrpc_transcript_counts_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "mcrpc_transcript_counts_metadata.csv")

# Create ExperimentHub metadata for Cao Esophageal transcript counts
cao_esophageal_transcript_counts_metadata = data.frame(
  Title = "cao_esophageal_transcript_counts",
  Description = "Transcript counts for all samples in cao_esophageal_wgbs_hg19 
  quantified using Kallisto with Gencode version 38 transcript annotation",
  BiocVersion = "3.18",
  Genome = "hg38",
  SourceType = "FASTQ",
  SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149609",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = NA,
  DataProvider = "University of California San Francisco",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "data.frame",
  DispatchClass = "FilePath",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "record/8422703/files/cao_esophageal_transcript_counts.tsv.gz"
)
write.csv(cao_esophageal_transcript_counts_metadata, "../extdata/cao_esophageal_transcript_counts_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "cao_esophageal_transcript_counts_metadata.csv")

# Create ExperimentHub metadata for TARGET rhabdoid transcript counts
target_rhabdoid_transcript_counts_metadata = data.frame(
  Title = "target_rhabdoid_transcript_counts",
  Description = "Transcript counts for 65 samples in target_rhabdoid_wgbs_hg19 
  quantified using Kallisto with Gencode version 38 transcript annotation",
  BiocVersion = "3.18",
  Genome = "hg38",
  SourceType = "FASTQ",
  SourceUrl = "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000470.v19.p8",
  SourceVersion = "1",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = NA,
  DataProvider = "Canada's Michael Smith Genome Sciences Centre",
  Maintainer = "Richard Heery <richardheery@gmail.com>",
  RDataClass = "data.frame",
  DispatchClass = "FilePath",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = "records/10019178/files/target_rhabdoid_kallisto_counts_all_transcripts.tsv.gz"
)
write.csv(target_rhabdoid_transcript_counts_metadata, "../extdata/target_rhabdoid_transcript_counts_metadata.csv", row.names = F)
ExperimentHubData::makeExperimentHubMetadata(pathToPackage = "~/git_repos/TumourMethData/", fileName = "target_rhabdoid_transcript_counts_metadata.csv")
