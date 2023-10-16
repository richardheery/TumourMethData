### Make cpgea_wgbs_hg38

# Load required packages
library(methodical)
library(dplyr)

# Get all CpG sites in hg38 as a GRanges
hg38_cpgs = methodical::extract_meth_sites_from_genome("BSgenome.Hsapiens.UCSC.hg38")

# Get paths to all bedGraph files
cpgea_bedgraphs = list.files("cpgea_bedgraphs", full.names = T)

# Create CPFEA sample annotation
cpgea_sample_metadata = data.frame(
  condition = ifelse(startsWith(basename(cpgea_bedgraphs), "T"), "Tumour", "Normal"),
  patient = as.character(readr::parse_number(basename(cpgea_bedgraphs))),
  row.names = gsub("_WGBS.CG.bg.gz", "", basename(cpgea_bedgraphs))
)

# Create methylation RSE from all bedGraphs. 
system.time({cpgea_meth_rse = methodical::make_meth_rse_from_bedgraphs(
  bedgraphs = cpgea_bedgraphs, meth_sites = hg38_cpgs, 
  sample_metadata = cpgea_sample_metadata, hdf5_dir = "cpgea_wgbs_hg38", 
  temporary_dir = "temp_hdf5_chunks", ncores = 20)})


### Make mcrpc_wgbs_hg38

# Get all downloaded methylation files
mcrpc_original_files = list.files("mcrpc_box_methylation_files", full.names = T)

# Get sample names from files
sample_names = gsub(".merged.CX_report.txt.gz", "", basename(mcrpc_original_files))
output_files = paste0("mcrpc_cpg_files/", sample_names, ".tsv.gz")

# Make a cluster
cl = makeCluster(10)
registerDoParallel(cl, 10)

# Filter original files for CpG sites using awk, convert chromosome names to UCSC format using sed and save to mcrpc_cpg_files directory. Took 15 minutes for one file. Took 2 hours for all files using 10 cores. 
system.time({foreach(file_number = seq_along(mcrpc_original_files))  %dopar% {
  system2("zcat", args = c(mcrpc_original_files[file_number], "|", "awk 'BEGIN {OFS=\"\t\"}; $6 == \"CG\"'", "|", "sed -f refseq_to_ucsc_chromosome_sed_commands.txt", "|", "gzip", ">", output_files[file_number]))
}})

# One of the original methylation files was named DTB-098-PRO2.merged.CX_report.txt.gz and so the cpg_file is renamed to DTB-098-PRO.tsv.gz
file.rename("mcrpc_cpg_files/DTB-098-PRO2.tsv.gz", "mcrpc_cpg_files/DTB-098-PRO.tsv.gz")

# Get paths to all CpG methylation files
mcrpc_cpg_files = list.files("mcrpc_cpg_files", full.names = T)

# Add sample names to filepaths
names(mcrpc_cpg_files) = gsub(".tsv.gz", "", basename(mcrpc_cpg_files))

# Download clinical data for MCRPC
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-020-0648-8/MediaObjects/41588_2020_648_MOESM3_ESM.xlsx", destfile = "supplementary_tables_zhao.xlsx")

# Load sample metadata. Samples can be mCRPC (73 samples), CMP (CpG methylator phenotype = 22 samples) or tSCNC (treatment-emergent small-cell neuroendocrine cancer = 5 samples)
mcrpc_sample_metadata = openxlsx::read.xlsx("supplementary_tables_zhao.xlsx", sheet = 1, startRow = 2)

# Change "mCRPC (CMP)" to just "CMP" in Type
mcrpc_sample_metadata$Type[mcrpc_sample_metadata$Type == "mCRPC (CMP)"] = "CMP"

# Select desired columns and save
mcrpc_sample_metadata = dplyr::select(mcrpc_sample_metadata, patient_id = Sample, 
    metastatis_site = Site, subtype = Type, age = Patient.Age)
mcrpc_sample_metadata$sex = "male"

# Remove BL and PRO from row.names
row.names(mcrpc_sample_metadata) = gsub("-", "_", gsub("-BL|-PRO", "", row.names(mcrpc_sample_metadata)))

# Change "-" to "_" in sample names
row.names(mcrpc_sample_metadata) = gsub("-", "_", row.names(mcrpc_sample_metadata))

# Put mcrpc_cpg_files in order of mcrpc_sample_metadata
mcrpc_cpg_files = mcrpc_cpg_files[row.names(mcrpc_sample_metadata)]

# Get reference CpG sites for hg38
hg38_cpgs_methrix = methrix::extract_CPGs("BSgenome.Hsapiens.UCSC.hg38")

# Create methrix object from CpG methylation files. Took 3 hours. 
system.time({mcrpc_methrix = methrix::read_bedgraphs(mcrpc_cpg_files, ref_cpgs = hg38_cpgs_methrix, 
  coldata = mcrpc_sample_metadata, chr_idx = 1, start_idx = 2, strand_idx = 3, M_idx = 4, U_idx = 5, 
  zero_based = F, stranded = T, collapse_strands = T, h5 = T, h5_dir = "mcrpc_methrix_h5")})

# Convert methrix into a meth RSE and resave
mcrpc_wgbs_hg38 = methodical::methrix_to_rse(mcrpc_methrix)
HDF5Array::quickResaveHDF5SummarizedExperiment(mcrpc_rse)
file.rename("mcrpc_methrix_h5", "mcrpc_wgbs_hg38/")

# Subset mcrpc_wgbs_hg38 for chromosome 11
mcrpc_wgbs_hg38_chr11 = mcrpc_wgbs_hg38[seqnames(mcrpc_wgbs_hg38) == "chr11", ]
HDF5Array::saveHDF5SummarizedExperiment(mcrpc_wgbs_hg38_chr11, "mcrpc_wgbs_hg38_chr11")

### Make tcga_wgbs_hg38

# Get paths to all TCGA WGBS tcga_bedgraphs and put in the same order as mcrpc_sample_metadata
tcga_bedgraphs = list.files("tcga_tcga_bedgraphs", full.names = T)
names(tcga_bedgraphs) = gsub("\\..*", "", basename(tcga_bedgraphs))

# Download metadata for TCGA samples
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0073-4/MediaObjects/41588_2018_73_MOESM3_ESM.xlsx", "tcga_metadata.xlsx")

# Get table with sample metadata. Seems there are two mistakes in the barcode column in the metadata, 
# with TCGA_COAD_3158 and TCGA_COAD_N3158 instead of TCGA_COAD_3518 and TCGA_COAD_N3518. This is corrected
tcga_sample_metadata = xlsx::read.xlsx("tcga_metadata.xlsx", sheetIndex = 1)
tcga_sample_metadata$barcode[tcga_sample_metadata$barcode == "TCGA_COAD_3158"] = "TCGA_COAD_3518"
tcga_sample_metadata$barcode[tcga_sample_metadata$barcode == "TCGA_COAD_N3158"] = "TCGA_COAD_N3518"

# Convert caseID to real TCGA barcodes
tcga_sample_metadata$real_barcodes = sapply(tcga_sample_metadata$caseID, function(x) 
  TCGAutils::UUIDtoBarcode(x, from_type = c("case_id"))$submitter_id)

# Replace fileName column by name of tcga_bedgraphs
tcga_tcga_sample_metadata$fileName = basename(tcga_bedgraphs[tcga_tcga_sample_metadata$barcode])

# Recreate sample metadata with submitter, sample_type, project and filename
tcga_sample_metadata = transmute(tcga_sample_metadata, submitter = real_barcodes, 
  sample_type = ifelse(tcga_sample_metadata$tumor.normal == "tumor", "01", "11"), project = cancer, filename = fileName)

# Add a column combining submitter and sample type 
tcga_sample_metadata$submitter_and_sample_type = paste(tcga_sample_metadata$submitter, tcga_sample_metadata$sample_type, sep = "_")

# Change "-" in submitter_and_sample_type to "_"
tcga_sample_metadata$submitter_and_sample_type = gsub("-", "_", tcga_sample_metadata$submitter_and_sample_type)

# Change submitter_and_sample_type to row.names
tcga_sample_metadata = tibble::column_to_rownames(tcga_sample_metadata, "submitter_and_sample_type")

# Put tcga_bedgraphs in the same order as tcga_sample_metadata
tcga_bedgraphs = tcga_bedgraphs[gsub("\\..*", "", tcga_sample_metadata$filename)]

# Create a HDF5-backed RangedSummarizedExperiment object for TCGA WGBS data. Took 40 minutes with 16 cores. 
system.time({tcga_wgbs_rse = methodical::make_meth_rse_from_tcga_bedgraphs(
  tcga_bedgraphs = tcga_bedgraphs, meth_sites = hg38_cpgs, tcga_sample_metadata = tcga_sample_metadata, 
  hdf5_dir = "tcga_wgbs_hg38", temporary_dir = "temp_hdf5_chunks", ncores = 16)})

### Make cao_esophageal_wgbs_hg19

# Get names of ESCA bedgraph files and extract sample names from them
cao_esophageal_bedgraphs = list.files("cao_esophageal_bedgraphs", full.names = T)
names(cao_esophageal_bedgraphs) = stringr::str_match(cao_esophageal_bedgraphs, "cao_esophageal_bedgraphs/G.*_(.*?)_CpG_methylation.bed.gz")[, 2]

# Load metadata 
cao_esophageal_sample_metadata = read.table("esca_cao_sample_metadata.tsv", header = T, row.names = 1)

# Ensure bedGraphs are in the same order as sample_metadata
cao_esophageal_bedgraphs = cao_esophageal_bedgraphs[row.names(sample_metadata)]

# Get hg19 CpGs
hg19_cpgs = methodical::extract_meth_sites_from_genome(genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)

# Create methrix object from CpG methylation files. Took 30 minutes with 5 cores. 
system.time({cao_esophageal_wgbs_hg19 = methodical::make_meth_rse_from_bedgraphs(bedgraphs = cao_esophageal_bedgraphs, convert_percentages = T, 
    meth_sites = hg19_cpgs, sample_metadata = cao_esophageal_sample_metadata, hdf5_dir = "cao_esophageal_wgbs_hg19", 
    temporary_dir = "temp_chunks", ncores = 5)})