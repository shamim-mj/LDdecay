#' Calculate Global Linkage Disequilibrium Decay and Render Regression Curves
#'
#' Evaluates genome-wide linkage disequilibrium (LD) decay profiles across populations. 
#' The function automates the structural synchronization of genotype matrices and variant map arrays,
#' filters out unanchored genomic scaffolds, processes pairwise r2 coefficients using local 
#' vectorized window blocks, trims the physical distance space to isolate decay plateaus, fits a 
#' non-linear least squares (NLS) regression model, and plots a profile.
#'
#' @param geno A data frame or matrix of genotype dosages. Must contain a `"taxa"` column 
#'   identifying individual samples, with all remaining columns containing variant IDs (scaled 0/1/2).
#' @param map A data frame containing variant coordinate details. Must include columns named
#'   `"SNP"`, `"Chromosome"`, and `"Position"`.
#' @param max_distance_bp Numeric. The maximum physical distance cutoff (in base-pairs) to trim the 
#'   long x-axis tail and focus the visualization window on where the curve plateaus. Default is \code{500000} (500kb).
#' @param max_plot_points Numeric. The maximum number of pairwise data points to sample and display 
#'   on the final scatterplot scatter matrix to prevent system RAM memory failure. Default is \code{100000}.
#' @param r2_threshold Numeric. The base-line horizontal intercept boundary showing target LD background limits. 
#'   Default is \code{0.2}.
#' @param output_image_path Character string pointing to a local file location (including file extension) 
#'   to save the rendered chart asset. Set to \code{NULL} to bypass saving. Default is \code{"Figures/ld_short_range.png"}.
#'
#' @return A list structure containing two key package elements:
#'   \item{plot}{The generated composite \code{ggplot} visualization object layout configuration.}
#'   \item{ld_data}{A filtered \code{data.frame} holding pairwise variant physical distances (\code{d}) and linkage values (\code{r2}), excluding self-correlations.}
#' @export
#'
#' @import ggplot2
#' @importFrom dplyr filter mutate select sample_n .data bind_rows tibble
#' @importFrom stringr str_detect
#' @importFrom scales comma_format
#' @importFrom stats cor
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
    dplyr::sample_n(plot_samples_count)
  
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
