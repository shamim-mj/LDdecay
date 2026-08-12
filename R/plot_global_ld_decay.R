#' Calculate and Plot Global Linkage Disequilibrium Decay
#'
#' Calculates genome-wide pairwise linkage disequilibrium (LD) as squared
#' Pearson correlations (\eqn{r^2}) between genotype markers within a
#' user-defined physical distance. The function synchronizes genotype and
#' marker-map data, removes unanchored scaffold markers, calculates pairwise
#' \eqn{r^2} values within chromosomes, and generates a scatterplot of LD
#' against physical distance. A nonlinear least-squares (NLS) regression
#' curve is fitted using the model \eqn{r^2 = 1/(1 + C d)}, where \eqn{d}
#' is physical distance and \eqn{C} is the fitted decay parameter.
#'
#' To reduce memory and rendering requirements, a maximum number of pairwise
#' observations can be randomly sampled for visualization. The returned
#' \code{ld_data} object contains the complete set of calculated pairwise
#' observations retained after filtering, whereas the plot is generated from
#' the sampled observations.
#'
#' @param geno A data frame or matrix containing genotype dosages. A column
#'   named \code{"taxa"} is required to identify samples, while the remaining
#'   columns should contain variant identifiers matching the \code{SNP}
#'   column in \code{map}. Genotype values should be coded numerically
#'   (e.g., 0, 1, and 2).
#'
#' @param map A data frame containing physical marker information. It must
#'   contain the columns \code{"SNP"}, \code{"Chromosome"}, and
#'   \code{"Position"}, corresponding to variant identifier, chromosome, and
#'   physical position in base pairs, respectively.
#'
#' @param max_distance_bp Numeric. Maximum physical distance in base pairs
#'   between two markers for which pairwise LD is calculated. Default is
#'   \code{500000} (500 kb).
#'
#' @param max_plot_points Numeric. Maximum number of pairwise LD observations
#'   randomly sampled for visualization. This limits the number of points
#'   rendered in the plot and helps reduce memory and plotting requirements.
#'   Default is \code{100000}.
#'
#' @param r2_threshold Numeric. Value used to draw a horizontal reference line
#'   on the plot, representing a user-defined LD threshold. Default is
#'   \code{0.2}.
#'
#' @param output_image_path Character string specifying the path, including
#'   file extension, where the generated plot should be saved. If \code{NULL},
#'   the plot is not saved to a file. Default is
#'   \code{"Figures/ld_short_range.png"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{A \code{ggplot} object showing pairwise \eqn{r^2} against
#'   physical distance, with an NLS LD-decay curve and user-defined
#'   \eqn{r^2} reference line.}
#'   \item{ld_data}{A data frame containing the calculated pairwise LD
#'   observations after removal of missing values and self-comparisons.
#'   Columns include \code{d}, the physical distance between marker pairs,
#'   and \code{r2}, the corresponding squared correlation coefficient.}
#' }
#'
#' @details
#' Markers classified as unanchored scaffolds are removed before LD
#' calculations. Only variants present in both the genotype and map datasets
#' are retained. Pairwise LD is calculated only between markers located on
#' the same chromosome and within \code{max_distance_bp} base pairs.
#'
#' The fitted decay model is:
#'
#' \deqn{
#' r^2 = \frac{1}{1 + C d}
#' }
#'
#' where \eqn{d} is physical distance and \eqn{C} is estimated using nonlinear
#' least squares. The regression curve is fitted to a random sample of at
#' most \code{max_plot_points} pairwise observations. The complete filtered
#' pairwise dataset is returned in \code{ld_data}.
#'
#' @examples
#' \dontrun{
#' library(LDdecay)
#'
#' # Calculate genome-wide LD decay
#' ld_result <- plot_global_ld_decay(
#'   geno = geno,
#'   map = map,
#'   max_distance_bp = 500000,
#'   max_plot_points = 100000,
#'   r2_threshold = 0.2
#' )
#'
#' # Display the LD-decay plot
#' ld_result$plot
#'
#' # Inspect the calculated pairwise LD data
#' head(ld_result$ld_data)
#'
#' # Save the plot to a custom location
#' ld_result <- plot_global_ld_decay(
#'   geno = geno,
#'   map = map,
#'   output_image_path = "Figures/global_ld_decay.png"
#' )
#' }
#'
#' @export
#'

#' 
plot_global_ld_decay <- function(geno, map, max_distance_bp = 500000, max_plot_points = 100000, 
                                 r2_threshold = 0.2, output_image_path = "Figures/ld_short_range.png") {
  
  old_scipen <- options(scipen = 999)
  on.exit(options(old_scipen), add = TRUE)
  
  message("Step 1: Synchronizing and cleaning genotype and map data frames...")
  
  if (!"taxa" %in% colnames(geno)) {
    stop("Error: The genotype dataset must contain a designated 'taxa' tracking column identifying samples.")
  }
  
  map_clean <- map %>%
    dplyr::filter(!stringr::str_detect(tolower(.data$Chromosome), "^scaffold"))
  
  common_snps <- intersect(colnames(geno)[colnames(geno) != "taxa"], map_clean$SNP)
  
  if (length(common_snps) < 2) {
    stop("Error: Insufficient variant overlap matching intersection bounds found between map arrays and genotype arrays.")
  }
  
  map_clean <- map_clean[map_clean$SNP %in% common_snps, ]
  map_clean <- map_clean[order(map_clean$Chromosome, map_clean$Position), ]
  common_snps_sorted <- map_clean$SNP
  
  message(paste(" -> Synchronized datasets: Processing", length(common_snps_sorted), "variants across", nrow(geno), "taxa samples."))
  
  geno_subset <- geno[, common_snps_sorted, drop = FALSE]
  geno_mat <- matrix(
    suppressWarnings(as.numeric(as.matrix(geno_subset))), 
    nrow = nrow(geno_subset), 
    ncol = ncol(geno_subset)
  )
  colnames(geno_mat) <- common_snps_sorted
  rownames(geno_mat) <- geno$taxa
  
  message("Step 2: Calculating sliding physical window pairwise r2 metrics...")
  ld_list <- list()
  n_markers <- length(common_snps_sorted)
  
  for (i in 1:(n_markers - 1)) {
    focal_pos <- map_clean$Position[i]
    focal_chr <- map_clean$Chromosome[i]
    
    j <- i + 1
    neighbors <- c()
    
    while (j <= n_markers && 
           map_clean$Chromosome[j] == focal_chr && 
           (map_clean$Position[j] - focal_pos) <= max_distance_bp) {
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
  
  message("Step 3: Post-processing calculation results...")
  ld_filtered <- dplyr::bind_rows(ld_list) %>% 
    dplyr::filter(!is.na(.data$r2) & .data$d != 0)
  
  if (nrow(ld_filtered) == 0) {
    stop("Error: No pairwise variant sets discovered within the specified physical distance window parameters.")
  }
  
  message("Step 4: Compiling visual elements and sampling data matrix limits...")
  plot_samples_count <- min(max_plot_points, nrow(ld_filtered))
  ld_sampled <- ld_filtered %>% 
    dplyr::slice_sample(n = plot_samples_count)
  
  message("Step 5: Fitting and plotting native NLS regression layer curve...")
  ld_decay_plot <- ggplot2::ggplot(ld_sampled, ggplot2::aes(x = .data$d, y = .data$r2)) +
    ggplot2::geom_point(alpha = 0.2, color = "gray50") +
    ggplot2::geom_smooth(
      method = "nls",
      formula = y ~ 1 / (1 + C * x),
      method.args = list(start = list(C = 0.00001)),
      color = "blue",
      linewidth = 1.2,
      se = FALSE
    ) +
    ggplot2::geom_hline(yintercept = r2_threshold, color = "black", linetype = "dashed") +
    ggplot2::scale_x_continuous(labels = scales::comma_format()) +
    ggplot2::labs(
      title = "Genome-Wide Global Linkage Disequilibrium (LD) Decay Profile",
      subtitle = paste("Model: r2 ~ 1/(1 + C*d) | Short-Range Max:", format(max_distance_bp, big.mark = ","), "bp"),
      x = "Physical Distance (bp)", 
      y = expression(r^{2}~value)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(family = "serif", size = 12, color = "black", face = "bold"),
      axis.title = ggplot2::element_text(family = "serif", size = 14, color = "black", face = "bold"),
      axis.text = ggplot2::element_text(family = "serif", size = 12, color = "black", face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
  
  if (!is.null(output_image_path)) {
    dir_name <- dirname(output_image_path)
    if (!dir.exists(dir_name) && dir_name != ".") {
      dir.create(dir_name, recursive = TRUE)
    }
    message(paste("-> Saving high-resolution publication plot asset to:", output_image_path))
    ggplot2::ggsave(
      filename = output_image_path, 
      plot = ld_decay_plot, 
      dpi = 300, 
      width = 8, 
      height = 4
    )
  }
  
  message("Done! Global analysis completed successfully.")
  
  return(list(
    plot = ld_decay_plot,
    ld_data = ld_filtered
  ))
}
