#' Download a WGBS methylation dataset from TumourMethData
#'
#' The HDF5 file and RDS file to construct a RangedSummarizedExperiment for the specified dataset
#' are downloaded into the ExperimentHub cache located at `ExperimentHub::getExperimentHubOption("cache")`.
#' A copy of the RDS file pointing to the HDF5 file in the cache is saved in the specified directory. 
#'
#' @param dataset Name of the dataset to download WGBS data from. Must be one of the datsets listed in data(TumourMethDatasets).
#' @param dir Directory in which to save a the copy of the RDS file pointing to the
#' location of the HDF5 file in the ExperimentHub cache. Default is the current directory.
#' @return A RangedSummarizedExperiment with methylation values from the specified dataset.
#' @export
#' @examples
#' mcrpc_wgbs_hg38_chr11 = TumourMethData::download_meth_dataset(dataset = "mcrpc_wgbs_hg38_chr11")
#' print(mcrpc_wgbs_hg38_chr11)
download_meth_dataset = function(dataset, dir = getwd()){
 
  # Load TumourMethDatasets
  data("TumourMethDatasets", package = "TumourMethData")
 
  # Check that dataset is one of the allowed options
  if(!dataset %in% TumourMethDatasets$dataset_name){
    stop("dataset should be one of the dataset names in TumourMethDatasets")
  }
 
  # Create output_dir from dir and dataset name
  output_file = paste0(dir, "/", dataset, ".rds")
 
  # Check if output_dir already exists
  if(!dir.exists(dir)){stop("dir doesn't exist")}
 
  # Extract the appropriate EH ID for the dataset
  eh_id = .experimenthub_ids[dataset, "wgbs"]
 
  # Create a connection to ExperimentHub and find the entry for the specified dataset
  eh  = ExperimentHub::ExperimentHub()
  dataset_files = eh[[eh_id]]
 
  # Check that two files were downloaded
  if(length(dataset_files) != 2){
    stop(paste("There were", length(dataset_files), "files downloaded, however there should be only 2 files per dataset"))
  }
 
  # Identify H5 file
  h5_file = dataset_files[which(sapply(dataset_files, function(x)
    tryCatch({rhdf5::h5ls(x); TRUE}, error = function(e) FALSE)))]
  rds_file = dataset_files[which(sapply(dataset_files, function(x)
    tryCatch({readRDS(x); TRUE}, error = function(e) FALSE)))]
 
  # Get the locatin of the ExperimentHubCache and change into the directory
  eh_cache = ExperimentHub::getExperimentHubOption("CACHE") 
  setwd(eh_cache)
  
  # Load rds_file_copy and change path of each assay to point to H5 file in ExperimentHub cache
  rse_temp = readRDS(rds_file)
  for(i in seq_along(SummarizedExperiment::assays(rse_temp, withDimnames = F))){
    slot(assay(rse_temp, i, withDimnames = F), "seed")@filepath = h5_file
  }
  saveRDS(rse_temp, output_file)
 
  # Create RangedSummarizedExperiment from output_dir
  rse = readRDS(output_file)
  return(rse)
 
}

#' Download a RNA-Seq counts dataset from TumourMethData
#'
#' @param dataset Name of the dataset to download. Must be one of the datasets 
#' listed in data(TumourMethDatasets) where `transcript_counts_available` is TRUE. 
#' @return A data.frame with RNA-Seq counts calcualted using Kallisto. 
#' @export
#' @examples
download_rnaseq_dataset = function(dataset){
  
  # Load TumourMethDatasets and filter for datasets with RNA-Seq
  data("TumourMethDatasets", package = "TumourMethData")
  TumourMethDatasets = TumourMethDatasets[TumourMethDatasets$transcript_counts_available, ]
  
  # Check that dataset is one of the allowed options
  if(!dataset %in% TumourMethDatasets$dataset_name){
    stop("dataset should be one of the dataset names in TumourMethDatasets where transcript_counts_available is TRUE")
  }
  
  # Extract the appropriate EH ID for the dataset
  eh_id = .experimenthub_ids[dataset, "rnaseq"]
  
  # Create a connection to ExperimentHub and find the entry for the specified dataset
  eh  = ExperimentHub::ExperimentHub()
  rnaseq_data = eh[[eh_id]]
  return(rnaseq_data)
  
}