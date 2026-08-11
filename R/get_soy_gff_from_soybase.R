#' Download, Decompress, and Import SoyBase GFF3 Annotation Files
#'
#' Automatically connects to the authoritative SoyBase repository data collection, 
#' downloads the compressed Williams 82 genome annotation layer (.gff3.gz), and 
#' extracts it into a memory-optimized data frame ready for variant structural mapping.
#'
#' @param dest_dir Character string pointing to the local directory where the downloaded 
#'   GFF3 file should be stored. Default is \code{"data"}.
#' @param overwrite Logical. If \code{TRUE}, overwrites the file if it already exists 
#'   locally in the target folder. Default is \code{FALSE}.
#'
#' @return A standard \code{data.frame} containing parsed genomic features (such as chromosomes, 
#'   start/end positions, strands, and associated metadata attributes). This function return gff3 file only for soybean. Download it manually for your crop of interest
#' @export
#'
#' @importFrom rtracklayer import
#' @importFrom utils download.file
#' 
#' @examples
#' \dontrun{
#' library(LDdecay)
#' 
#' # Downloads from SoyBase, decompresses, and builds your matrix cleanly
#' gff_table <- get_soy_gff_from_soybase(dest_dir = "data")
#' 
#' # Preview the results instantly!
#' head(gff_table)
#' }
#' 

get_soy_gff_from_soybase <- function(dest_dir = "data",overwrite = FALSE) {
  
  # Ensure the target directory structure exists
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }
  
  # Define authoritative SoyBase target paths and file names
  base_url  <- "https://data.soybase.org/Glycine/max/annotations/Wm82.gnm1.ann1.DvBy/"
  file_name <- "glyma.Wm82.gnm1.ann1.DvBy.gene_models_main.gff3.gz"
  full_url  <- paste0(base_url, file_name)
  
  dest_file_gz  <- file.path(dest_dir, file_name)
  dest_file_gff <- file.path(dest_dir, gsub("\\.gz$", "", file_name))
  
  # Check if uncompressed file already exists to save time and bandwidth
  if (file.exists(dest_file_gff) && !overwrite) {
    message("-> Found uncompressed GFF3 file locally. Bypassing download pipeline.")
    return(LDdecay::get_gff(dest_file_gff))
  }
  
  # Download compressed archive
  message("-> Initiating network download connection from SoyBase...")
  message("Source URL: ", full_url)
  
  utils::download.file(
    url      = full_url, 
    destfile = dest_file_gz, 
    mode     = "wb", 
    quiet    = FALSE
  )
  
  if (!file.exists(dest_file_gz)) {
    stop("Error: File transfer failed. Check your internet connection or the server URL path status.")
  }
  
  message("-> File download complete. Decompressing archive layer...")
  # Decompress natively within R using gzfile connections
  gz_con  <- gzfile(dest_file_gz, "rb")
  out_con <- file(dest_file_gff, "wb")
  
  writeBin(readBin(gz_con, "raw", n = 5e7), out_con)
  
  close(gz_con)
  close(out_con)
  
  # Remove the temporary downloaded compressed .gz archive file to save local storage space
  if (file.exists(dest_file_gz)) {
    file.remove(dest_file_gz)
  }
  
  message("-> Decompression successful! Saved clean file to: ", dest_file_gff)
  
  # Pass straight through to your existing package function to format into a light data frame
  return(LDdecay::get_gff(dest_file_gff))
}
