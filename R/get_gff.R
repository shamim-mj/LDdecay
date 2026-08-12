#' Import and Convert GFF3 Annotation Files to Standard Data Frames
#'
#' Reads a local genomic annotation file in GFF3 format and converts it into a
#' standard R data frame for downstream genomic analyses.
#'
#' @param gff_path Character string specifying the path to a GFF3 annotation file.
#'
#' @return A \code{data.frame} containing parsed genomic features, including
#'   chromosome, start and end coordinates, strand, feature type, and
#'   annotation attributes.
#'
#' @examples
#' \dontrun{
#' library(LDdecay)
#'
#' gff_table <- get_gff("data/annotation.gff3")
#' head(gff_table)
#' }
#'
#' @export
get_gff <- function(gff_path) {
  
  message("-> Loading GFF3 annotation layers into memory...")
  
  gff_data <- rtracklayer::import(
    gff_path,
    format = "gff3"
  )
  
  gff_data <- as.data.frame(gff_data)
  return(gff_data)
}
