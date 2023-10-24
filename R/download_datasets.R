#' Download a WGBS methylation dataset from TumourMethData
#' 
#' The HDF5 file and RDS file to construct a RangedSummarizedExperiment for the specified dataset 
#' are downloaded into the ExperimentHub cache located at `ExperimentHub::getExperimentHubOption("cache")`
#' and symbolic links to these are created in the specified directory by `dir`. 
#'
#' @param dataset Name of the dataset to download WGBS data from. Must be one of the datsets listed in data(TumourMethDatasets). 
#' @param dir Parent directory to create links to the HDF5 SummarizedExperiment dataset. 
#' A subdirectory with the dataset name will be created within this directory. Default is tempdir().
#' If a directory named `paste(dir, dataset, sep = "/")` already exists, 
#' will attempt to load a HDF5 summarized experiment from there and return an error if this fails.
#' @return A RangedSummarizedExperiment with methylation values from the specified dataset. 
#' @export
#' @examples
#' mcrpc_wgbs_hg38_chr11 = TumourMethData::download_meth_dataset(dataset = "mcrpc_wgbs_hg38_chr11")
#' print(mcrpc_wgbs_hg38_chr11)
download_meth_dataset = function(dataset, dir = tempdir()){
  
  # Load TumourMethDatasets
  data("TumourMethDatasets", package = "TumourMethData")
  
  # Check that dataset is one of the allowed options
  if(!dataset %in% TumourMethDatasets$dataset_name){
    stop("dataset should be one of the dataset names in TumourMethDatasets")
  }
  
  # Create output_dir from dir and dataset name
  output_dir = paste(dir, dataset, sep = "/")
  
  # Check if output_dir already exists
  if(!dir.exists(dir)){stop("dir doesn't exist")}
  if(dir.exists(output_dir)){
    tryCatch({
      rse = HDF5Array::loadHDF5SummarizedExperiment(output_dir)
      print(paste("A HDF5 SummarizedExperiment is already present in", output_dir, "and is being returned"))
      return(rse)
    }, 
      error = function(err) stop(paste(output_dir, "already exists but it does not contain a HDF5 SummarizedExperiment"))
  )}
  
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
  
  # Check that one H5 file and one RDS file were found
  if(length(h5_file) == 0){stop("No HDF5 file was downloaded")}
  if(length(rds_file) == 0){stop("No RDS file was downloaded")}
  
  # Create directory and create symlinks to files in it
  dir.create(output_dir)
  h5_link = paste(output_dir, "assays.h5", sep = "/")
  rds_link = paste(output_dir, "se.rds", sep = "/")
  R.utils::createLink(link = h5_link, target = h5_file)
  R.utils::createLink(link = rds_link, target = rds_file)
  
  # Create RangedSummarizedExperiment from output_dir
  rse = HDF5Array::loadHDF5SummarizedExperiment(output_dir)
  return(rse)
  
}

#' Download an RNA-Seq counts dataset from TumourMethData
#' 
#' A TSV file with the RNA-Seq counts is downloaded into the ExperimentHub cache 
#' located at `ExperimentHub::getExperimentHubOption("cache")` and is read into R as a data.frame. 
#'
#' @param dataset Name of the dataset to download. Must be one of the datasets 
#' listed in data(TumourMethDatasets) where `transcript_counts_available` is TRUE. 
#' @return A data.frame with RNA-Seq counts calculated using Kallisto. 
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
  rnaseq_counts_file = eh[[eh_id]]
  rnaseq_data = read.table(rnaseq_counts_file, sep = "\t", row.names = 1)
  return(rnaseq_data)
  
}
