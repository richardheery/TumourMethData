# Create a methylation RSE for all CPGEA WGBS bedGraphs

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
  hdf5_dir = "cpgea_meth_rse_hdf5", temporary_dir = "temp_hdf5_chunks", ncores = 20
  )})