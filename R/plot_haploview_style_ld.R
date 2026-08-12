#' Generate a Haploview-Style Inverted Triangle LD Heatmap
#'
#' Calculates pairwise linkage disequilibrium (LD) as squared Pearson
#' correlations (\eqn{r^2}) among a selected set of variants and renders
#' the resulting matrix as a Haploview-style LD heatmap. The heatmap can be
#' displayed in an inverted-triangle configuration and is aligned with the
#' physical positions of the selected variants.
#'
#' The function supports several ways of specifying the variants to display,
#' including explicit SNP identifiers, numeric row indices, and character
#' ranges such as \code{"1:60"}. If no subset is provided, the first 60
#' variants in the map are selected automatically.
#'
#' @param geno A data frame or matrix containing genotype dosages. Columns
#'   corresponding to the variants specified in \code{snp_subset} must be
#'   present in the genotype data. Genotypes should be coded numerically
#'   (e.g., 0, 1, and 2).
#'
#' @param map A data frame containing physical marker information. It must
#'   contain the columns \code{"SNP"} and \code{"Position"}, corresponding
#'   to variant identifiers and physical positions in base pairs,
#'   respectively. A \code{"Chromosome"} column may also be present for
#'   compatibility with standard marker-map files.
#'
#' @param snp_subset Specifies the variants to include in the heatmap.
#'   The argument can be:
#'   \itemize{
#'     \item a character vector of SNP identifiers;
#'     \item a numeric vector of row indices in \code{map}; or
#'     \item a single character range such as \code{"1:60"}.
#'   }
#'   If \code{NULL}, the first 60 variants in \code{map} are selected.
#'
#' @param title Character string specifying the title displayed on the
#'   heatmap. Default is \code{"Pairwise LD"}.
#'
#' @param color_palette Character vector defining the colors used to
#'   represent LD values. The default is \code{grDevices::heat.colors(20)}.
#'
#' @param flip Logical. If \code{TRUE}, the LDheatmap is displayed in the
#'   inverted-triangle orientation. Default is \code{TRUE}.
#'
#' @return An object returned by \code{LDheatmap::LDheatmap}, containing the
#'   rendered LD heatmap and associated grid graphics components.
#'
#' @details
#' The function first identifies the requested variants from the marker map
#' and synchronizes them with the genotype matrix. Pairwise LD is calculated
#' as:
#'
#' \deqn{
#' r^2 = cor(G_i, G_j)^2
#' }
#'
#' where \eqn{G_i} and \eqn{G_j} are genotype dosage vectors for two
#' variants. Pairwise correlations are calculated using
#' \code{use = "pairwise.complete.obs"}, allowing available non-missing
#' genotype observations to be used for each pair.
#'
#' The resulting square \eqn{r^2} matrix is supplied to
#' \code{LDheatmap::LDheatmap}, together with the physical positions and
#' variant identifiers. The function therefore provides a convenient way
#' to visualize local LD structure for a selected group of variants.
#'
#' @examples
#' \dontrun{
#' library(LDdecay)
#'
#' # Plot the first 60 variants in the marker map
#' ld_plot <- plot_haploview_style_ld(
#'   geno = geno,
#'   map = map
#' )
#'
#' # Plot variants 1 through 30
#' ld_plot <- plot_haploview_style_ld(
#'   geno = geno,
#'   map = map,
#'   snp_subset = "1:30"
#' )
#'
#' # Plot a specific set of SNPs
#' selected_snps <- c(
#'   "ss715620779",
#'   "ss715620801",
#'   "ss715620845"
#' )
#'
#' ld_plot <- plot_haploview_style_ld(
#'   geno = geno,
#'   map = map,
#'   snp_subset = selected_snps,
#'   title = "Candidate SNP LD"
#' )
#'
#' # Display the heatmap
#' ld_plot
#' }
#'
#' @export
#'
#' @importFrom stats cor

plot_haploview_style_ld <- function(geno, map, snp_subset = NULL, title = "Pairwise LD", 
                              color_palette = grDevices::heat.colors(20), flip = TRUE) {
  
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
