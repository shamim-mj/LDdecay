#' Plot Linkage Disequilibrium Heatmap Around a Focal Variant
#'
#' Isolates flanking markers inside a physical window, filters out variants with 
#' a low Minor Allele Frequency (MAF), and renders a correlation matrix heatmap.
#'
#' @param geno A numeric matrix or data frame of genotype dosages. Rows are samples, columns are SNPs.
#' @param map A data frame containing variant metadata map records.
#' @param focal_snp Character string identifying the central target SNP.
#' @param window_bp Numeric. Distance upstream/downstream to include variants. Default is \code{100000} (100kb).
#' @param marker_col Character. Column name for variant IDs in the map file. Default is \code{"SNP"}.
#' @param chr_col Character. Column name for chromosomes in the map file. Default is \code{"Chromosome"}.
#' @param pos_col Character. Column name for base-pair positions in the map file. Default is \code{"Position"}.
#' @param min_maf Numeric. Minor Allele Frequency threshold to filter variants. Default is \code{0.05}.
#'
#' @return A list containing the \code{ggplot} object, the raw $r^2$ correlation matrix, and filtered maps/genotypes.
#' @export
#'
#' @import ggplot2
#' @importFrom stats cor
#' @importFrom dplyr tibble .data
plot_snp_ld <- function(geno, map, focal_snp, window_bp = 100000, 
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

