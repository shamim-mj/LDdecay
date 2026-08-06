#' Import and Convert GFF3 Annotation Files to Standard Data Frames
#'
#' Reads a local genomic annotation file in GFF3 format and converts it into a 
#' standard R data frame table, making structural genomic features accessible for 
#' downstream data operations.
#'
#' @param gff_path A character string pointing to the local file path of a valid `.gff3` file.
#'
#' @return A standard `data.frame` containing parsed genomic features (such as chromosomes, 
#'   start/end positions, strands, and associated metadata attributes).
#' @export
#'
#' @importFrom rtracklayer import
get_gff <- function(gff_path) {
  
  message("-> Loading GFF3 annotation layers into memory...")
  gff_data <- rtracklayer::import(gff_path, format = "gff3")
  
  # Conversion to standard R data frame table
  gff_table <- as.data.frame(gff_data)
  
  return(gff_table)
}
