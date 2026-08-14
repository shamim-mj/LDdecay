#' Calculate and Plot Binned Empirical Linkage Disequilibrium Decay
#'
#' Calculates pairwise linkage disequilibrium (LD) as squared Pearson
#' correlation (\eqn{r^2}) between genotype markers and summarizes LD across
#' physical-distance intervals. The resulting empirical LD decay profile is
#' visualized as the mean \eqn{r^2} across user-defined distance bins.
#'
#' Genotype and physical-map data are first synchronized using the intersection
#' of marker identifiers present in both datasets. Unanchored scaffold
#' chromosomes are excluded, and markers are ordered by chromosome and physical
#' position. Pairwise \eqn{r^2} values are then calculated between markers on
#' the same chromosome within the specified maximum physical distance.
#'
#' @param geno A data frame or matrix containing genotype dosages. The data
#'   must contain a column named \code{"taxa"} identifying individual samples,
#'   with all remaining columns corresponding to variant identifiers. Genotype
#'   values should be numeric dosage values such as 0, 1, and 2.
#' @param map A data frame containing physical coordinates for the variants.
#'   The data frame must contain columns named \code{"SNP"},
#'   \code{"Chromosome"}, and \code{"Position"}, corresponding to variant
#'   identifiers, chromosome assignments, and physical positions in base pairs,
#'   respectively.
#' @param bin_size_bp Numeric value specifying the width of the physical-distance
#'   bins used to summarize pairwise LD. The value is expressed in base pairs.
#'   Default is \code{10000} (10 kb).
#' @param max_distance_bp Numeric value specifying the maximum physical distance
#'   over which pairwise LD is calculated and displayed. The value is expressed
#'   in base pairs. Default is \code{500000} (500 kb).
#' @param r2_threshold Numeric value specifying the \eqn{r^2} threshold displayed
#'   as a horizontal dashed reference line on the plot. Default is \code{0.2}.
#' @param output_image_path Character string specifying the file path at which
#'   the generated plot should be saved. The file extension determines the
#'   output format supported by \code{\link[ggplot2]{ggsave}}. Set to
#'   \code{NULL} to return the plot without saving an image. Default is
#'   \code{"Figures/ld_binned_decay.png"}.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{plot}{A \code{ggplot} object containing the empirical LD decay
#'   profile.}
#'   \item{binned_data}{A data frame containing the summarized LD values for
#'   each physical-distance bin. Columns include \code{bin}, the lower
#'   boundary of the distance bin in base pairs; \code{mean_r2}, the mean
#'   pairwise \eqn{r^2}; \code{median_r2}, the median pairwise \eqn{r^2}; and
#'   \code{count}, the number of pairwise comparisons contributing to the bin.}
#' }
#'
#' @details
#' The function calculates LD separately within each chromosome and does not
#' calculate LD between markers located on different chromosomes. Markers
#' located on chromosomes whose names begin with \code{"scaffold"} are
#' excluded before analysis.
#'
#' Only markers present in both \code{geno} and \code{map} are retained.
#' Markers are ordered by chromosome and physical position before pairwise
#' comparisons are performed.
#'
#' For each marker, the function identifies downstream markers on the same
#' chromosome within \code{max_distance_bp}. Pairwise LD is calculated as the
#' squared Pearson correlation between genotype dosage vectors using
#' pairwise-complete observations.
#'
#' Pairwise LD values are grouped into physical-distance bins of width
#' \code{bin_size_bp}. For each bin, the function calculates the mean
#' \eqn{r^2}, median \eqn{r^2}, and number of contributing marker pairs.
#' The resulting mean \eqn{r^2} values are used to construct the empirical
#' LD decay curve.
#'
#' A small calculation buffer is added internally to the requested maximum
#' distance to reduce edge effects when assigning pairwise observations to
#' the final distance bins.
#'
#' If \code{output_image_path} is not \code{NULL}, the output directory is
#' created automatically when necessary and the plot is saved at 300 dpi.
#'
#' @examples
#' \dontrun{
#' library(LDdecay)
#'
#' # Calculate empirical LD decay using 10-kb distance bins
#' ld_decay <- plot_binned_ld_global_decay(
#'   geno = geno,
#'   map = map,
#'   bin_size_bp = 10000,
#'   max_distance_bp = 500000
#' )
#'
#' # Display the LD decay plot
#' ld_decay$plot
#'
#' # Inspect the summarized LD data
#' head(ld_decay$binned_data)
#'
#' # Run the analysis without saving an image
#' ld_decay <- plot_binned_ld_global_decay(
#'   geno = geno,
#'   map = map,
#'   output_image_path = NULL
#' )
#' }
#'
#' @export
#'
#' @import ggplot2
#' @importFrom dplyr filter mutate group_by summarise select .data
#' @importFrom stringr str_detect
#' @importFrom scales comma_format
#' @importFrom stats cor median
#' @importFrom dplyr filter mutate group_by summarise select .data
#' @importFrom magrittr  %>%

plot_binned_ld_global_decay <- function(geno,
                                        map,
                                        bin_size_bp = 10000,
                                        max_distance_bp = 500000, 
                                        r2_threshold = 0.2,
                                        output_image_path = "Figures/ld_binned_decay.png") {
  
  # Enforce non-scientific text layouts across function boundaries
  old_scipen <- options(scipen = 999)
  on.exit(options(old_scipen), add = TRUE)
  
  message("Step 1: Synchronizing and cleaning genotype and map data frames...")
  
  if (!"taxa" %in% colnames(geno)) {
    stop("Error: The genotype dataset must contain a designated 'taxa' tracking column identifying samples.")
  }
  
  # Drop structural unanchored scaffolds from physical map structures
  map_clean <- map %>%
    dplyr::filter(!stringr::str_detect(tolower(.data$Chromosome), "^scaffold"))
  
  # Isolate overlapping variant intersection markers present across both datasets
  common_snps <- intersect(colnames(geno)[colnames(geno) != "taxa"], map_clean$SNP)
  
  if (length(common_snps) < 2) {
    stop("Error: Insufficient variant overlap matching intersection bounds found between map arrays and genotype arrays.")
  }
  
  # Synchronize frames and explicitly sort by Chromosome and Position for window checks
  map_clean <- map_clean[map_clean$SNP %in% common_snps, ]
  map_clean <- map_clean[order(map_clean$Chromosome, map_clean$Position), ]
  common_snps_sorted <- map_clean$SNP
  
  message(paste(" -> Synchronized datasets: Processing", length(common_snps_sorted), "variants across", nrow(geno), "taxa samples."))
  
  # Convert genotypes to standard calculation matrix layout
  geno_mat <- as.matrix(geno[, common_snps_sorted, drop = FALSE])
  storage.mode(geno_mat) <- "numeric"
  
  # Add a calculation buffer zone to avoid edge clipping effects on the final bin
  calculation_buffer_bp <- bin_size_bp * 2
  max_calculation_window <- max_distance_bp + calculation_buffer_bp
  
  message(paste("Step 2: Calculating sliding physical window pairwise r2 metrics up to", max_calculation_window, "bp..."))
  ld_list <- list()
  n_markers <- length(common_snps_sorted)
  
  for (i in 1:(n_markers - 1)) {
    focal_pos <- map_clean$Position[i]
    focal_chr <- map_clean$Chromosome[i]
    
    j <- i + 1
    neighbors <- c()
    
    # Fast window check matching chromosome boundaries and maximum buffered physical distance thresholds
    while (j <= n_markers && 
           map_clean$Chromosome[j] == focal_chr && 
           (map_clean$Position[j] - focal_pos) <= max_calculation_window) {
      neighbors <- c(neighbors, j)
      j <- j + 1
    }
    
    if (length(neighbors) > 0) {
      r2_vals <- stats::cor(geno_mat[, i], geno_mat[, neighbors], use = "pairwise.complete.obs")^2
      
      ld_list[[i]] <- dplyr::tibble(
        d  = map_clean$Position[neighbors] - focal_pos,
        r2 = as.vector(r2_vals)
      )
    }
  }
  
  message("Step 3: Compiling and grouping distances into intervals...")
  ld_raw <- dplyr::bind_rows(ld_list) %>% 
    dplyr::filter(!is.na(.data$r2) & .data$d != 0)
  
  if (nrow(ld_raw) == 0) {
    stop("Error: No pairwise variant sets discovered within the specified physical distance window parameters.")
  }
  
  message(paste("-> Aggregating pairwise correlations into", bin_size_bp / 1000, "kb grouping intervals..."))
  
  # Assign markers to user-defined distance bins and aggregate summary statistics
  ld_binned <- ld_raw %>%
    dplyr::mutate(bin = floor(.data$d / bin_size_bp) * bin_size_bp) %>%
    # Filter the data layer down BEFORE grouping to isolate data from boundary artifacts
    dplyr::filter(.data$bin <= max_distance_bp) %>%
    dplyr::group_by(.data$bin) %>%
    dplyr::summarise(
      mean_r2   = mean(.data$r2, na.rm = TRUE),
      median_r2 = stats::median(.data$r2, na.rm = TRUE),
      count     = dplyr::n(),
      .groups   = "drop"
    )
  
  # Step 4: Build plot structure configurations using full package namespaces
  message("-> Rendering empirical summary visualization profile charts...")
  binned_plot <- ggplot2::ggplot(ld_binned, ggplot2::aes(x = .data$bin, y = .data$mean_r2)) +
    ggplot2::geom_line(color = "red", linewidth = 1) +
    ggplot2::geom_point(color = "darkred") +
    # Switch out scale limits for coord_cartesian to avoid line clipping bugs entirely
    ggplot2::coord_cartesian(xlim = c(0, max_distance_bp)) +
    ggplot2::scale_x_continuous(
      n.breaks = 10, 
      labels   = scales::comma_format()
    ) +
    ggplot2::geom_hline(yintercept = r2_threshold, color = "black", linetype = "dashed") +
    ggplot2::labs(
      title    = "Empirical Binned Linkage Disequilibrium Decay Profile",
      subtitle = paste("Grouping Size Radius:", format(bin_size_bp, big.mark = ","), "bp Intervals"),
      x        = "Physical Distance (bp)", 
      y        = expression(r^{2}~value)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text        = ggplot2::element_text(family = "serif", size = 12, color = "black", face = "bold"),
      axis.title  = ggplot2::element_text(family = "serif", size = 14, color = "black", face = "bold"),
      axis.text   = ggplot2::element_text(family = "serif", size = 8, color = "black", face = "bold"),
      plot.title  = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
  
  # Step 5: Handle image file generation pathways safely
  if (!is.null(output_image_path)) {
    dir_name <- dirname(output_image_path)
    if (!dir.exists(dir_name) && dir_name != ".") {
      dir.create(dir_name, recursive = TRUE)
    }
    message(paste("-> Saving publication-ready summary line chart asset to:", output_image_path))
    ggplot2::ggsave(
      filename = output_image_path, 
      plot     = binned_plot, 
      dpi      = 300, 
      width    = 8, 
      height   = 4
    )
  }
  
  message("Done! Binned analysis complete.")
  return(list(plot = binned_plot, binned_data = ld_binned))
}
