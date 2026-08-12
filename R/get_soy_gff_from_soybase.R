#' Download the SoyBase Williams 82 GFF3 Annotation
#'
#' Downloads the Williams 82 soybean genome annotation from SoyBase,
#' decompresses the GFF3 archive, and converts the annotation to a standard
#' data frame using \code{\link{get_gff}}.
#'
#' This function is specifically designed for soybean (\emph{Glycine max})
#' and retrieves the Williams 82 genome annotation from SoyBase.
#'
#' If the uncompressed GFF3 file already exists in \code{dest_dir} and
#' \code{overwrite = FALSE}, the existing file is used and no download is
#' performed.
#'
#' @param dest_dir Character string specifying the directory in which the
#'   GFF3 annotation file will be downloaded and stored. The directory is
#'   created recursively if it does not already exist. Default is
#'   \code{"data"}.
#' @param overwrite Logical indicating whether an existing uncompressed GFF3
#'   file should be replaced. Default is \code{FALSE}.
#'
#' @return A \code{data.frame} containing the parsed Williams 82 soybean
#'   genome annotation.
#'
#' @details
#' The function downloads
#' \code{glyma.Wm82.gnm1.ann1.DvBy.gene_models_main.gff3.gz}
#' from SoyBase, decompresses the archive, removes the compressed file,
#' and passes the resulting GFF3 file to \code{\link{get_gff}} for parsing.
#'
#' @references
#' SoyBase. Soybean Genome Database.
#' \url{https://www.soybase.org/}
#'
#' @examples
#' \dontrun{
#' library(LDdecay)
#'
#' gff_table <- get_soy_gff_from_soybase(
#'   dest_dir = "data"
#' )
#'
#' head(gff_table)
#' }
#'
#' @export
#' @importFrom utils download.file
get_soy_gff_from_soybase <- function(
    dest_dir = "data",
    overwrite = FALSE
) {
  
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }
  
  base_url <- paste0(
    "https://data.soybase.org/",
    "Glycine/max/annotations/",
    "Wm82.gnm1.ann1.DvBy/"
  )
  
  file_name <- "glyma.Wm82.gnm1.ann1.DvBy.gene_models_main.gff3.gz"
  
  full_url <- paste0(base_url, file_name)
  
  dest_file_gz <- file.path(
    dest_dir,
    file_name
  )
  
  dest_file_gff <- file.path(
    dest_dir,
    sub("\\.gz$", "", file_name)
  )
  
  if (file.exists(dest_file_gff) && !overwrite) {
    
    message(
      "-> Found uncompressed GFF3 file locally. ",
      "Skipping download."
    )
    
    return(
      get_gff(dest_file_gff)
    )
  }
  
  message("-> Initiating download from SoyBase...")
  message("Source URL: ", full_url)
  
  utils::download.file(
    url = full_url,
    destfile = dest_file_gz,
    mode = "wb"
  )
  
  if (!file.exists(dest_file_gz)) {
    stop(
      "File download failed. Check your internet connection ",
      "or the SoyBase URL."
    )
  }
  
  message(
    "-> File download complete. Decompressing archive..."
  )
  
  R.utils::gunzip(
    filename = dest_file_gz,
    destname = dest_file_gff,
    overwrite = overwrite,
    remove = FALSE
  )
  
  if (file.exists(dest_file_gz)) {
    file.remove(dest_file_gz)
  }
  
  message(
    "-> Decompression successful! Saved clean file to: ",
    dest_file_gff
  )
  
  get_gff(dest_file_gff)
}