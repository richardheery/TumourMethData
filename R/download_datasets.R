#' Download one of the datasets from TumourMethData
#'
#' @param dataset Name of the dataset to download. Must be one of the datsets listed in data(TumourMethDatasets). 
#' @param dir Parent directory to save the HDF5 SummarizedExperiment datasets. A directory with the dataset name 
#' will be created within this directory and populated with the RDS and HDF5 files constituting the dataset.  
#' @return A RangedSummarizedExperiment with methylation values from the specified dataset. 
#' @export
#' @examples
#' tcga_wgbs_hg38 = TumourMethData::download_dataset(dataset = "tcga_wgbs_hg38", dir = ".")
#' print(tcga_wgbs_hg38)
download_dataset = function(dataset, dir = "."){
  
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
  if(dir.exists(output_dir)){stop(
    paste(output_dir, "already exists")
  )}
  
  # Create a connection to ExperimentHub and find the entry for the specified dataset
  eh  = ExperimentHub::ExperimentHub()
  dataset_entry = BiocGenerics::subset(eh, title == dataset)
  dataset_files = dataset_entry[[1]]
  
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
  
  # Create directory and move files into it
  dir.create(paste(dir, dataset, sep = "/"))
  file.rename(h5_file, paste(output_dir, "assays.h5", sep = "/"))
  file.rename(rds_file, paste(output_dir, "se.rds", sep = "/"))
  
  # Create RangedSummarizedExperiment from output_dir
  rse = HDF5Array::loadHDF5SummarizedExperiment(output_dir)
  return(rse)
  
}