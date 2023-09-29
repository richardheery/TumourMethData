### Make cpgea_wgbs_hg38

# Load required packages
library(methodicalFinal)

# Set working directory to CPGEA WGBS directory on scratch
setwd("~/mounts/hpcnfs_mount/scratch/MS/cpgea/wgbs")

# Get paths to all bedGraph files
bgs = list.files("all_bedgraphs_hg38", full.names = T)

# Create sample annotation
sample_annotation = data.frame(
  condition = ifelse(startsWith(basename(bgs), "T"), "Tumour", "Normal"),
  patient = as.character(readr::parse_number(basename(bgs))),
  row.names = gsub("_WGBS.CG.bg.gz", "", basename(bgs))
)

# Get all CpG sites in hg38 as a GRanges
hg38_cpgs = methodicalFinal:::cpg_genome_ranges_hg38

# Create methylation RSE from all bedGraphs. Took 5 hours. 
system.time({cpgea_meth_rse = methodicalFinal::make_meth_rse_from_bedgraphs(
  bedgraphs = bgs, meth_sites = hg38_cpgs, sample_metadata = sample_annotation, 
  hdf5_dir = "cpgea_wgbs_hg38", temporary_dir = "temp_hdf5_chunks", ncores = 20
  )})

### Make tcga_wgbs_hg38

# Get paths to all TCGA WGBS bedgraphs and put in the same order as sample_metadata
bedgraphs = list.files("tcga_bedgraphs", full.names = T)
names(bedgraphs) = gsub("\\..*", "", basename(bedgraphs))

# Replace fileName column by name of bedGraphs
sample_metadata$fileName = basename(bedgraphs[sample_metadata$barcode])

# Recreate sample metadata with submitter, sample_type, project and filename
sample_metadata = transmute(sample_metadata, submitter = real_barcodes, 
  sample_type = ifelse(sample_metadata$tumor.normal == "tumor", "01", "11"), project = cancer, filename = fileName)

# Add a column combining submitter and sample type 
sample_metadata$submitter_and_sample_type = paste(sample_metadata$submitter, sample_metadata$sample_type, sep = "_")

# Change "-" in submitter_and_sample_type to "_"
sample_metadata$submitter_and_sample_type = gsub("-", "_", sample_metadata$submitter_and_sample_type)
write.table(sample_metadata, "tcga_wgbs_sample_metadata.tsv", sep = "\t", row.names = F, quote = F)

# Change submitter_and_sample_type to row.names
sample_metadata = data.table::fread("tcga_wgbs_sample_metadata.tsv", sep = "\t")
sample_metadata = tibble::column_to_rownames(sample_metadata, "submitter_and_sample_type")

# Put bedgraphs in the same order as sample_metadata
bedgraphs = bedgraphs[gsub("\\..*", "", sample_metadata$filename)]

# Get all CpG sites in hg38 as a GRanges
hg38_cpgs = methodicalFinal:::cpg_genome_ranges_hg38

# Create a HDF5-backed RangedSummarizedExperiment object for TCGA WGBS data. Took 40 minutes with 16 cores. 
system.time({tcga_wgbs_rse = methodical::make_meth_rse_from_bedgraphs(
  bedgraphs = bedgraphs, meth_sites = hg38_cpgs, sample_metadata = sample_metadata, 
  hdf5_dir = "tcga_wgbs_hg38", temporary_dir = "temp_hdf5_chunks", ncores = 16)})