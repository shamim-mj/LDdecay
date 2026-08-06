#' Generate a Haploview-Style Inverted Triangle LD Heatmap
#'
#' Extracts a targeted subset of structural variants, calculates a square pairwise 
#' correlation matrix ($r^2$), and renders a publication-quality Haploview-style 
#' inverted triangle heatmap aligned with physical chromosome positions.
#'
#' @param geno A data frame or matrix of genotype dosages. Must contain a tracking column 
#'   (e.g., `"taxa"`) with all remaining columns containing variant column names (scaled 0/1/2).
#' @param map A data frame containing physical variant metadata map records. Must include 
#'   columns named `"SNP"` and `"Position"`.
#' @param snp_subset A vector specifying which SNPs to plot. Can be a character vector of SNP names, 
#'   a numeric vector of row indices (e.g., `1:60`), or a string range (e.g., `"1:60"`). 
#'   Default is `NULL` (automatically selects the first 60 variants).
#' @param title Character string defining the primary header title of the heatmap layout. Default is `"Pairwise LD"`.
#' @param color_palette A vector of colors specifying the heatmap grading scale. Default is `heat.colors(20)`.
#' @param flip Logical. If `TRUE`, displays the standard inverted triangle configuration format. Default is `TRUE`.
#'
#' @return An object layout of class \code{LDheatmap} containing structural grid graphics configurations.
#' @export
#'
#' @importFrom stats cor
#' @importFrom dplyr filter
#' @importFrom LDheatmap LDheatmap
plot_haploview_ld <- function(geno, map, snp_subset = NULL, title = "Pairwise LD", 
                              color_palette = heat.colors(20), flip = TRUE) {
  
  message("Step 1: Validating inputs and filtering variant coordinates...")
  
  # Ensure map is a standard data frame
  map <- as.data.frame(map)
  
  # Step 1b: Smart parsing of the snp_subset parameter
  if (is.null(snp_subset)) {
    message(" -> No subset provided. Automatically isolating the first 60 variants available in the map file.")
    num_snps <- min(60, nrow(map))
    target_snps <- as.character(map$SNP[1:num_snps])
    
  } else if (is.character(snp_subset) && length(snp_subset) == 1 && grepl("^\\d+:\\d+$", snp_subset)) {
    # Handles string ranges like "1:60" safely
    message(paste(" -> Detected string range expression:", snp_subset, ". Parsing into row indices."))
    bounds <- as.numeric(unlist(strsplit(snp_subset, ":")))
    row_indices <- bounds[1]:bounds[2]
    row_indices <- row_indices[row_indices >= 1 & row_indices <= nrow(map)]
    target_snps <- as.character(map$SNP[row_indices])
    
  } else if (is.numeric(snp_subset)) {
    # Handles numeric ranges like 1:60 directly
    message(" -> Detected numeric row indices.")
    valid_indices <- snp_subset[snp_subset >= 1 & snp_subset <= nrow(map)]
    target_snps <- as.character(map$SNP[valid_indices])
    
  } else {
    # Handles standard list vectors of explicit SNP names
    message(" -> Detected character vector of explicit variant identifiers.")
    target_snps <- as.character(snp_subset)
  }
  
  # Filter map to include only the target SNPs and maintain order integrity
  small_map <- map[map$SNP %in% target_snps, , drop = FALSE]
  small_map <- small_map[match(target_snps, small_map$SNP), ]
  small_map <- small_map[!is.na(small_map$SNP), ]
  
  final_snps <- as.character(small_map$SNP)
  
  if (length(final_snps) < 2) {
    stop("Error: Insufficient target variants (fewer than 2) found across your dataset coordinates. Check your row indices or SNP names.")
  }
  
  # Ensure target SNPs exist as columns in the genotype data
  final_snps <- intersect(final_snps, colnames(geno))
  small_map <- small_map[small_map$SNP %in% final_snps, ]
  
  if (length(final_snps) < 2) {
    stop("Error: Fewer than 2 matching variants were found within the genotype column names.")
  }
  
  message(paste("Step 2: Subsetting genotypes for", length(final_snps), "synchronized variants..."))
  
  # Extract numeric dosages
  raw_subset <- geno[, final_snps, drop = FALSE]
  raw_subset <- as.matrix(raw_subset)
  storage.mode(raw_subset) <- "numeric"
  
  message("Step 3: Calculating square pairwise r2 correlation matrix...")
  ld_square_matrix <- stats::cor(raw_subset, use = "pairwise.complete.obs")^2
  
  message("Step 4: Executing LDheatmap layout render loops...")
  
  # Call the LDheatmap engine
  heatmap_output <- LDheatmap::LDheatmap(
    gdat               = ld_square_matrix,
    genetic.distances  = small_map$Position,
    title              = title,
    SNP.name           = small_map$SNP,
    color              = color_palette,
    flip               = flip
  )
  
  message("Done! Triangle map rendered successfully.")
  return(invisible(heatmap_output))
}
