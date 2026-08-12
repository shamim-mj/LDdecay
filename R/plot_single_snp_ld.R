#' Plot Linkage Disequilibrium Heatmap Around a Focal Variant
#'
#' Calculates pairwise linkage disequilibrium (LD) as squared Pearson
#' correlations (\eqn{r^2}) among variants surrounding a focal SNP and
#' visualizes the resulting LD matrix as a heatmap. Variants are restricted
#' to the same chromosome as the focal SNP, filtered by a user-defined
#' physical window, and screened using a minimum minor allele frequency (MAF)
#' threshold.
#'
#' The focal SNP is highlighted with an asterisk in the heatmap axis labels.
#' The resulting plot provides a visual representation of local LD structure
#' around a significant or otherwise prioritized variant.
#'
#' @param geno A data frame or matrix containing genotype dosages, with rows
#'   representing samples and columns representing variants. If present, a
#'   column named \code{"taxa"} is automatically removed before LD
#'   calculations. Genotype values should be numerically coded, such as
#'   0, 1, and 2.
#'
#' @param map A data frame containing physical marker information. It must
#'   contain columns identifying the variant, chromosome, and physical
#'   position. The corresponding column names are specified using
#'   \code{marker_col}, \code{chr_col}, and \code{pos_col}.
#'
#' @param focal_snp Character string specifying the focal variant around
#'   which the LD region is constructed. The focal SNP must be present in
#'   both the marker map and genotype data.
#'
#' @param window_bp Numeric. Physical distance in base pairs upstream and
#'   downstream of the focal SNP to include in the analysis. Default is
#'   \code{100000} (100 kb on each side of the focal SNP).
#'
#' @param marker_col Character string specifying the column in \code{map}
#'   containing variant identifiers. Default is \code{"SNP"}.
#'
#' @param chr_col Character string specifying the column in \code{map}
#'   containing chromosome identifiers. Default is \code{"Chromosome"}.
#'
#' @param pos_col Character string specifying the column in \code{map}
#'   containing physical marker positions in base pairs. Default is
#'   \code{"Position"}.
#'
#' @param min_maf Numeric. Minimum minor allele frequency required for a
#'   variant to be retained in the LD analysis. Variants with MAF below this
#'   threshold are removed. Default is \code{0.05}.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{A \code{ggplot} object showing the pairwise local LD
#'   matrix as a heatmap.}
#'   \item{r2}{A square matrix containing pairwise squared Pearson
#'   correlation coefficients (\eqn{r^2}) among the retained variants.}
#'   \item{map}{A data frame containing the marker metadata for variants
#'   retained after regional and MAF filtering.}
#'   \item{genotypes}{A numeric matrix containing genotype dosages for the
#'   retained variants.}
#' }
#'
#' @details
#' The function first identifies the chromosome and physical position of the
#' focal SNP. It then selects markers located on the same chromosome and
#' within \code{window_bp} base pairs upstream or downstream of the focal
#' variant.
#'
#' Variants not present in the genotype matrix are removed. Minor allele
#' frequency is calculated from genotype dosage values as:
#'
#' \deqn{
#' AF = \frac{\mathrm{mean}(G)}{2}
#' }
#'
#' and minor allele frequency is calculated as:
#'
#' \deqn{
#' MAF = \min(AF, 1 - AF).
#' }
#'
#' Variants with MAF below \code{min_maf} are excluded before calculating
#' pairwise LD.
#'
#' Pairwise LD is calculated as the squared Pearson correlation between
#' genotype dosage vectors:
#'
#' \deqn{
#' r^2 = cor(G_i, G_j)^2.
#' }
#'
#' Correlations are calculated using \code{use = "pairwise.complete.obs"},
#' allowing each pairwise comparison to use available non-missing genotype
#' observations.
#'
#' The heatmap displays the lower triangular portion of the pairwise
#' correlation matrix. The focal SNP is identified by an asterisk in the
#' axis labels.
#'
#' @examples
#' \dontrun{
#' library(LDdecay)
#'
#' # Plot LD surrounding a focal SNP using a 100-kb window
#' ld_result <- plot_single_snp_ld(
#'   geno = geno,
#'   map = map,
#'   focal_snp = "ss715620779",
#'   window_bp = 100000
#' )
#'
#' # Display the heatmap
#' ld_result$plot
#'
#' # Inspect the pairwise LD matrix
#' head(ld_result$r2)
#'
#' # Inspect the markers retained in the analysis
#' head(ld_result$map)
#' }
#'
#' @export
#'
#' 
plot_single_snp_ld <- function(geno, map, focal_snp, window_bp = 100000, 
                        marker_col = "SNP", chr_col = "Chromosome", 
                        pos_col = "Position", min_maf = 0.05) {
  
  # Cleanly drop the character 'taxa' tracking column before processing numbers
  if ("taxa" %in% colnames(geno)) {
    geno <- geno[, colnames(geno) != "taxa", drop = FALSE]
  }
  
  # Safe matrix extraction template to completely bypass coercion warnings
  geno_mat <- matrix(
    suppressWarnings(as.numeric(as.matrix(geno))), 
    nrow = nrow(geno), 
    ncol = ncol(geno)
  )
  colnames(geno_mat) <- colnames(geno)
  
  if (!focal_snp %in% map[[marker_col]]) stop("Focal variant missing from map data records.")
  if (!focal_snp %in% colnames(geno_mat)) stop("Focal variant column missing from genotype matrix.")
  
  # Extract focal metrics
  focal_row <- map[map[[marker_col]] == focal_snp, , drop = FALSE]
  focal_chr <- focal_row[[chr_col]][1]
  focal_pos <- focal_row[[pos_col]][1]
  
  # Subset and sort regional map records
  region_map <- map[
    map[[chr_col]] == focal_chr & 
      map[[pos_col]] >= (focal_pos - window_bp) & 
      map[[pos_col]] <= (focal_pos + window_bp), , drop = FALSE
  ]
  
  region_map <- region_map[region_map[[marker_col]] %in% colnames(geno_mat), , drop = FALSE]
  region_map <- region_map[order(region_map[[pos_col]]), , drop = FALSE]
  markers    <- as.character(region_map[[marker_col]])
  
  if (length(markers) < 2) stop("Fewer than two operational markers found within requested window constraints.")
  
  # Filter genotypes by MAF values
  region_geno <- geno_mat[, markers, drop = FALSE]
  af  <- colMeans(region_geno, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  
  keep        <- is.finite(maf) & maf >= min_maf
  region_geno <- region_geno[, keep, drop = FALSE]
  region_map  <- region_map[keep, , drop = FALSE]
  markers     <- colnames(region_geno)
  
  if (length(markers) < 2) stop("Insufficient variant count remaining following rare variant MAF processing loops.")
  
  # Pairwise correlation matrix calculations
  r2 <- stats::cor(region_geno, use = "pairwise.complete.obs")^2
  
  # Pivot to long format matrix quickly
  plot_data <- as.data.frame(as.table(r2))
  names(plot_data) <- c("Marker1", "Marker2", "r2")
  plot_data$x <- match(plot_data$Marker1, markers)
  plot_data$y <- match(plot_data$Marker2, markers)
  plot_data   <- plot_data[plot_data$y <= plot_data$x, , drop = FALSE]
  
  # Apply functional asterisk labels to the focal SNP
  marker_labels <- markers
  marker_labels[marker_labels == focal_snp] <- paste0(focal_snp, " *")
  
  # Construct ggplot container assets
  ld_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$r2)) +
    ggplot2::geom_tile(color = "grey80", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(
      low = "white", high = "red3", limits = c(0, 1), 
      na.value = "grey90", name = expression(r^2)
    ) +
    ggplot2::scale_x_continuous(breaks = seq_along(markers), labels = marker_labels, position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_along(markers), labels = marker_labels, expand = c(0, 0)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste("LD Matrix Map Flanking", focal_snp),
      subtitle = paste0("Chromosome ", focal_chr, ": ", format(focal_pos - window_bp, big.mark = ","), " to ", format(focal_pos + window_bp, big.mark = ","), " bp"),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5, size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(face = "bold")
    )
  
  return(list(plot = ld_plot, r2 = r2, map = region_map, genotypes = region_geno))
}

