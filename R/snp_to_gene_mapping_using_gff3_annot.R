#' Map SNPs to Nearby Genomic Features using GFF3 Annotations
#'
#' Scans a physical window upstream and downstream of target SNPs to find overlapping
#' gene models extracted from a standard GFF3 file using high-performance non-equi joins.
#'
#' @param snp_list A character vector of target SNP identifiers to query. e.g. list of all significant SNPs detected in GWAS
#' @param gff_path A character string pointing to the local paths of a standard \code{.gff3} file. Downloaded from https://data.soybase.org/Glycine/max/annotations/Wm82.gnm1.ann1.DvBy/
#' @param map_data A data frame containing physical map coordinates. Must include columns such as
#'   \code{SNP}, \code{Chromosome}, and \code{Position}.
#' @param window_bp Numeric. The search radius in base-pairs flanking the variant. Default is \code{50000} (50kb).
#'
#' @return A polished \code{tibble} containing mapped variants, distance to midpoints, orientation, and descriptions.
#' @export
#'
#' @importFrom rtracklayer import
#' @importFrom dplyr mutate filter select arrange tibble case_when .data inner_join
#' @importFrom stringr str_extract
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' library(LDdecay)
#' library(tidyverse)
#' # Define target test SNPs
#' test_snps <- c("ss715620779", "ss715620778", "ss715620770")

#' # Run mapping pipeline
#' test_mapping <- snp_to_gene_mapping_using_gff3_annot(
#' snp_list = test_snps,
#' gff_path = path/to/gff3/file, 
#' map_data = map, 
#' window_bp = 50000)

#' # View performance tables
#' print(test_mapping)
#' }
#' 
snp_to_gene_mapping_using_gff3_annot <- function(snp_list, gff_path, map_data, window_bp = 50000) {
  
  # Step 1: Clean and filter map data to target SNPs only
  snp_list <- trimws(as.character(snp_list))
  
  map_clean <- map_data %>%
    dplyr::mutate(
      SNP = trimws(as.character(.data$SNP)),
      # Force chromosome names to standard short format (e.g. "Gm03")
      Chromosome = stringr::str_extract(trimws(as.character(.data$Chromosome)), "Gm\\d+"),
      Position = as.numeric(.data$Position)
    ) %>%
    dplyr::filter(.data$SNP %in% snp_list)
  
  if (nrow(map_clean) == 0) {
    stop("Error: Zero target SNPs matched your physical map file.")
  }
  
  message("-> Found ", nrow(map_clean), " out of ", length(snp_list), " requested SNPs in the map file.")
  
  # Step 2: Import GFF3 and instantly strip down metadata to avoid memory inflation
  message("-> Importing GFF3 feature records...")
  gff_granges <- rtracklayer::import(gff_path, format = "gff3")
  # Conversion to standard R data frame table
  gff_granges <- as.data.frame(gff_granges)
  gff_genes <- gff_granges[gff_granges$type == "gene", ]
  
  if (length(gff_genes) == 0) {
    stop("Error: No structural feature tags labeled 'gene' found in the GFF3 file.")
  }
  
  message("-> Stripping and converting structural elements to minimal tables...")
  
  # Safe extraction of attributes without calling un-vectorized df inflation
  gff_meta <- S4Vectors::mcols(gff_genes)
  
  desc_vector <- if ("Note" %in% colnames(gff_meta)) {
    sapply(gff_meta$Note, function(x) {
      if (is.null(x) || length(x) == 0) "" else paste(unlist(x), collapse = "; ")
    })
  } else {
    ""
  }
  
  # Build a bare, ultra-light data frame for matching
  gff_df <- dplyr::tibble(
    chromosome  = stringr::str_extract(as.character(GenomicRanges::seqnames(gff_genes)), "Gm\\d+"),
    seqnames    = as.character(GenomicRanges::seqnames(gff_genes)),
    start       = as.numeric(GenomicRanges::start(gff_genes)),
    end         = as.numeric(GenomicRanges::end(gff_genes)),
    strand      = as.character(GenomicRanges::strand(gff_genes)),
    type        = "gene", # Safely hardcoded since the subset contains only genes
    gene_id     = if ("ID" %in% colnames(gff_meta)) as.character(gff_meta$ID) else NA_character_,
    gene_symbol = if ("Name" %in% colnames(gff_meta)) as.character(gff_meta$Name) else NA_character_,
    description = desc_vector
  ) %>% 
    dplyr::filter(!is.na(.data$chromosome)) # Drop scaffold fragments missing clear 'Gm' marks
  
  # Step 3: Run safe positional mapping via standard tables
  message("-> Calculating regional windows and matching locations...")
  
  final_table <- map_clean %>%
    # Initial quick join purely on matching chromosome ids
    dplyr::inner_join(gff_df, by = c("Chromosome" = "chromosome"), relationship = "many-to-many") %>%
    # Vectorized boundary filtering matching your exact logic
    dplyr::filter(.data$start <= (.data$Position + window_bp) & .data$end >= (.data$Position - window_bp)) %>%
    dplyr::mutate(
      gene_midpoint = (.data$start + .data$end) / 2,
      distance_to_snp = round(.data$gene_midpoint - .data$Position),
      orientation   = dplyr::case_when(
        .data$distance_to_snp > 0 ~ "Downstream",
        .data$distance_to_snp < 0 ~ "Upstream",
        TRUE ~ "Overlapping"
      )
    ) %>%
    dplyr::select(
      query_snp = .data$SNP,
      snp_chr = .data$Chromosome,
      snp_pos = .data$Position,
      .data$gene_id,
      .data$gene_symbol,
      .data$seqnames,
      .data$start,
      .data$end,
      .data$strand,
      .data$type,
      .data$distance_to_snp,
      .data$orientation,
      .data$description
    ) %>%
    dplyr::arrange(.data$query_snp, abs(.data$distance_to_snp))
  
  message("Done! Master table generated with ", nrow(final_table), " total entries.")
  return(final_table)
}
