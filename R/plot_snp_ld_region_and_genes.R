#' Plot Linkage Disequilibrium and Gene Annotations for a Genomic Region
#'
#' Generates a publication-quality composite visualization of linkage
#' disequilibrium (LD) and gene annotations surrounding a focal variant.
#' The function calculates SNP-specific pairwise LD around the focal SNP,
#' extracts gene models from a GFF3-derived annotation table, identifies
#' candidate genes using a user-defined regular expression, and combines
#' the results into a two-panel genomic region plot.
#'
#' The upper panel displays gene models positioned according to their
#' physical genomic coordinates. Genes matching \code{gene_pattern} are
#' highlighted and can be labeled automatically. Significant SNPs can
#' also be displayed as vertical reference lines with short labels.
#' The lower panel displays pairwise LD as an inverted triangular
#' representation, where the x-axis corresponds to physical genomic
#' position and the y-axis represents half the physical distance between
#' marker pairs.
#'
#' This function is designed for investigating candidate genomic regions
#' surrounding significant GWAS markers and can be used to evaluate
#' whether significant variants are located within or near genes with
#' potentially relevant functional annotations.
#'
#' @param map A data frame containing physical marker information. It must
#'   contain columns identifying the marker, chromosome, and physical
#'   position. The corresponding column names are specified by
#'   \code{marker_col}, \code{chr_col}, and \code{pos_col}.
#'
#' @param geno A data frame or matrix containing genotype dosages.
#'   Rows represent individual samples and columns represent genetic
#'   markers. A \code{"taxa"} column identifying individuals may be
#'   included and is automatically removed before LD calculations.
#'   Genotype values should be numeric dosages, typically coded as
#'   0, 1, and 2.
#'
#' @param gff_table A data frame containing genomic annotations parsed
#'   from a GFF3 file. The table must contain at least the columns
#'   \code{seqnames}, \code{start}, \code{end}, \code{strand}, and
#'   \code{type}. Gene annotations should be available as records with
#'   \code{type == "gene"}.
#'
#' @param chromosome Character string or numeric value identifying the
#'   chromosome to analyze. Chromosome identifiers are standardized to
#'   soybean-style names such as \code{"Gm03"}.
#'
#' @param focal_snp Character string identifying the focal or significant
#'   SNP around which SNP-specific LD is calculated. The SNP must be
#'   present in both \code{map} and \code{geno}.
#'
#' @param ld_window Numeric. Physical distance in base pairs upstream and
#'   downstream of \code{focal_snp} used when calculating SNP-specific LD.
#'   Default is \code{360000} (360 kb).
#'
#' @param region_start Numeric. Lower physical coordinate defining the
#'   genomic region displayed in the final plot. If \code{NULL}, the
#'   minimum position among available markers in the selected chromosome
#'   is used.
#'
#' @param region_end Numeric. Upper physical coordinate defining the
#'   genomic region displayed in the final plot. If \code{NULL}, the
#'   maximum position among available markers in the selected chromosome
#'   is used.
#'
#' @param marker_col Character string specifying the column in \code{map}
#'   containing unique marker identifiers. Default is \code{"SNP"}.
#'
#' @param chr_col Character string specifying the column in \code{map}
#'   containing chromosome identifiers. Default is \code{"Chromosome"}.
#'
#' @param pos_col Character string specifying the column in \code{map}
#'   containing physical base-pair positions. Default is \code{"Position"}.
#'
#' @param significant_snps Character vector containing SNP identifiers
#'   that should be highlighted in the genomic region plot. Significant
#'   SNPs are displayed as vertical dashed lines and assigned short labels
#'   (e.g., \code{S1}, \code{S2}). Default is \code{NULL}.
#'
#' @param gene_pattern Character string containing a regular expression
#'   used to identify candidate genes from their available annotation
#'   fields, including gene names, descriptions, notes, and ontology
#'   information. Matching is case-insensitive. The default pattern
#'   searches for ABC/PDR transporter-related annotations.
#'
#' @param gene_name Character string defining the functional category name
#'   used to label genes matching \code{gene_pattern} in the plot legend.
#'   Default is \code{"ABC/PDR Transport"}.
#'
#' @param label_candidate_genes Logical. If \code{TRUE}, labels genes matching
#'   \code{gene_pattern} above their corresponding gene models. Despite
#'   the historical argument name, this option is not restricted to
#'   cytochrome P450 genes and can be used with any gene annotation
#'   pattern. Default is \code{TRUE}.
#'
#' @param max_ld_distance Numeric. Optional maximum physical distance in
#'   base pairs between marker pairs retained for the LD panel. If
#'   \code{NULL}, all pairwise comparisons within the selected region are
#'   retained. Default is \code{NULL}.
#'
#' @param label_char_width Numeric. Scaling factor used to estimate the
#'   horizontal space occupied by gene labels when assigning labels to
#'   non-overlapping lanes. Larger values allocate more horizontal space
#'   to each label. Default is \code{0.010}.
#'
#' @param label_lane_spacing Numeric. Vertical spacing between successive
#'   gene-label lanes used to reduce label overlap. Default is \code{0.55}.
#'
#' @return A named list containing the composite plot and supporting
#'   genomic data:
#'   \describe{
#'     \item{\code{plot}}{Complete composite \code{patchwork} plot containing
#'       the gene annotation and LD panels.}
#'     \item{\code{gene_plot}}{The gene annotation panel as a
#'       \code{ggplot} object.}
#'     \item{\code{ld_plot}}{The LD triangular panel as a
#'       \code{ggplot} object.}
#'     \item{\code{genes}}{All gene models overlapping the selected genomic
#'       region, including gene classes, tracks, and label positions.}
#'     \item{\code{candidate_genes}}{Subset of \code{genes} matching
#'       \code{gene_pattern}. The element name is retained for backward
#'       compatibility and is not restricted to CYP genes.}
#'     \item{\code{marker_map}}{Marker map records retained within the
#'       selected genomic region.}
#'     \item{\code{significant_map}}{Map records corresponding to
#'       \code{significant_snps}.}
#'     \item{\code{r2}}{Square pairwise LD correlation matrix for markers
#'       retained in the selected region.}
#'     \item{\code{ld_data}}{Long-format pairwise LD data containing marker
#'       identities, physical positions, pairwise \eqn{r^2}, marker
#'       distance, and plotting coordinates.}
#'     \item{\code{chromosome}}{Standardized chromosome identifier used
#'       for the analysis.}
#'     \item{\code{region_start}}{Lower physical coordinate of the plotted
#'       region.}
#'     \item{\code{region_end}}{Upper physical coordinate of the plotted
#'       region.}
#'   }
#'
#' @details
#' The function performs the following workflow:
#'
#' \enumerate{
#'   \item Calculates SNP-specific pairwise LD around the focal SNP using
#'         \code{\link{plot_single_snp_ld}}.
#'   \item Validates the marker map and GFF3-derived annotation table.
#'   \item Standardizes chromosome identifiers to the \code{Gm01},
#'         \code{Gm02}, ..., \code{Gm20} format.
#'   \item Selects markers located on the requested chromosome and within
#'         the specified genomic region.
#'   \item Converts the pairwise LD matrix into long format and retains
#'         one triangular half of the pairwise comparisons.
#'   \item Extracts gene models overlapping the selected genomic region.
#'   \item Searches gene annotation fields using \code{gene_pattern} to
#'         identify candidate genes.
#'   \item Assigns genes to non-overlapping horizontal tracks.
#'   \item Assigns candidate-gene labels to non-overlapping vertical lanes.
#'   \item Adds user-defined significant SNPs as reference lines and labels.
#'   \item Constructs a shared physical-position scale between the gene and
#'         LD panels.
#'   \item Combines the gene and LD panels using \code{patchwork}.
#' }
#'
#' Candidate genes are identified using a regular expression applied to
#' available annotation fields such as \code{Name}, \code{ID},
#' \code{Note}, \code{Description}, \code{Dbxref}, and
#' \code{Ontology_term}. Therefore, users can adapt the function to
#' investigate different gene families or biological pathways by changing
#' \code{gene_pattern} and \code{gene_name}.
#'
#' The LD panel uses the physical midpoint of each marker pair as its
#' horizontal coordinate and half of the physical marker distance as its
#' vertical coordinate. This produces a triangular representation similar
#' to commonly used LD-region visualizations.
#'
#' @examples
#' \dontrun{
#' library(LDdecay)
#'
#' # Define the focal genomic region
#' focal_snp <- "ss715620779"
#' chromosome <- "Gm15"
#'
#' # Define candidate-gene annotation pattern
#' gene_pattern <- "ABC-2/plant PDR ABC transporter|ABC|PDR transporter"
#' gene_name <- "ABC/PDR Transport"
#'
#' # Generate LD and gene annotation plot
#' result <- plot_snp_ld_region_and_genes(
#'   map = data(map),
#'   geno = data(geno),
#'   gff_table = gff_table,
#'   chromosome = chromosome,
#'   focal_snp = focal_snp,
#'   ld_window = 360000,
#'   gene_pattern = gene_pattern,
#'   gene_name = gene_name
#' )
#'
#' # Display the composite plot
#' result$plot
#'
#' # View candidate genes
#' result$candidate_genes
#'
#' # View pairwise LD data
#' head(result$ld_data)
#' }
#'
#' @seealso
#' \code{\link{plot_single_snp_ld}},
#' \code{\link{get_gff}}
#'
#' @import ggplot2
#' @importFrom dplyr filter mutate select
#' @importFrom stringr str_detect
#' @importFrom rlang .data
#' @importFrom stats cor setNames
#' @importFrom scales label_number breaks_pretty
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom grid unit
#' @export

plot_snp_ld_region_and_genes <- function(
    map,
    geno,
    gff_table,
    chromosome, 
    focal_snp,
    ld_window,
    region_start = NULL,
    region_end = NULL, 
    marker_col = "SNP",
    chr_col = "Chromosome",
    pos_col = "Position", 
    significant_snps = NULL, 
    gene_pattern = "ABC-2/plant PDR ABC transporter|ABC|PDR transporter", 
    gene_name = "ABC/PDR Transport", 
    label_candidate_genes = TRUE, max_ld_distance = NULL, 
    label_char_width = 0.010, label_lane_spacing = 0.55
) {
  
  #--------------------------------------------------#
  # 1. Check required packages
  #--------------------------------------------------#
  
  required_packages <- c(
    "ggplot2",
    "patchwork",
    "scales"
  )
  
  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]
  
  if (length(missing_packages) > 0) {
    stop(
      "Install the following packages first: ",
      paste(missing_packages, collapse = ", ")
    )
  }
  
  #--------------------------------------------------#
  # calculate SNP specific LD
  #--------------------------------------------------#
  
  ld_result <- plot_single_snp_ld(
    geno = geno,
    map = map,
    focal_snp = focal_snp, # change for other snps, 
    window_bp = ld_window,
    min_maf = 0.05
  )
  r2 <- ld_result$r2
  
  #--------------------------------------------------#
  # 2. Validate inputs
  #--------------------------------------------------#
  required_map_columns <- c(
    marker_col,
    chr_col,
    pos_col
  )
  
  if (!all(required_map_columns %in% names(map))) {
    stop(
      "The map is missing one or more required columns: ",
      paste(required_map_columns, collapse = ", ")
    )
  }
  
  required_gff_columns <- c(
    "seqnames",
    "start",
    "end",
    "strand",
    "type"
  )
  
  if (!all(required_gff_columns %in% names(gff_table))) {
    stop(
      "The GFF table is missing one or more required columns: ",
      paste(required_gff_columns, collapse = ", ")
    )
  }
  
  if (is.null(rownames(r2)) || is.null(colnames(r2))) {
    stop(
      "The r2 matrix must have marker names as both row names ",
      "and column names."
    )
  }
  
  if (
    !is.numeric(label_char_width) ||
    length(label_char_width) != 1 ||
    label_char_width <= 0
  ) {
    stop("'label_char_width' must be one positive number.")
  }
  
  if (
    !is.numeric(label_lane_spacing) ||
    length(label_lane_spacing) != 1 ||
    label_lane_spacing <= 0
  ) {
    stop("'label_lane_spacing' must be one positive number.")
  }
  
  #--------------------------------------------------#
  # 3. Standardize chromosome names
  #--------------------------------------------------#
  
  standardize_chromosome <- function(x) {
    
    x <- as.character(x)
    
    chromosome_number <- gsub(
      "[^0-9]",
      "",
      x
    )
    
    chromosome_number <- suppressWarnings(
      as.integer(chromosome_number)
    )
    
    ifelse(
      is.na(chromosome_number),
      NA_character_,
      sprintf("Gm%02d", chromosome_number)
    )
  }
  
  chromosome_short <- standardize_chromosome(
    chromosome
  )[1]
  
  if (is.na(chromosome_short)) {
    stop(
      "Could not interpret chromosome: ",
      chromosome
    )
  }
  
  map <- as.data.frame(map)
  gff_table <- as.data.frame(gff_table)
  
  map$chr_short <- standardize_chromosome(
    map[[chr_col]]
  )
  
  gff_table$chr_short <- sub(
    "^.*\\.(Gm[0-9]+)$",
    "\\1",
    as.character(gff_table$seqnames)
  )
  
  gff_table$chr_short <- standardize_chromosome(
    gff_table$chr_short
  )
  
  #--------------------------------------------------#
  # 4. Prepare the marker map
  #--------------------------------------------------#
  
  map[[marker_col]] <- as.character(
    map[[marker_col]]
  )
  
  map[[pos_col]] <- as.numeric(
    map[[pos_col]]
  )
  
  available_markers <- intersect(
    rownames(r2),
    colnames(r2)
  )
  
  region_map <- map[
    map$chr_short == chromosome_short &
      map[[marker_col]] %in% available_markers,
    ,
    drop = FALSE
  ]
  
  if (nrow(region_map) < 2) {
    stop(
      "Fewer than two markers from the LD matrix were found on ",
      chromosome_short,
      "."
    )
  }
  
  if (is.null(region_start)) {
    region_start <- min(
      region_map[[pos_col]],
      na.rm = TRUE
    )
  }
  
  if (is.null(region_end)) {
    region_end <- max(
      region_map[[pos_col]],
      na.rm = TRUE
    )
  }
  
  if (
    !is.finite(region_start) ||
    !is.finite(region_end) ||
    region_start >= region_end
  ) {
    stop(
      "'region_start' must be smaller than 'region_end'."
    )
  }
  
  region_map <- region_map[
    region_map[[pos_col]] >= region_start &
      region_map[[pos_col]] <= region_end,
    ,
    drop = FALSE
  ]
  
  region_map <- region_map[
    order(region_map[[pos_col]]),
    ,
    drop = FALSE
  ]
  
  # Remove duplicated marker records if present.
  region_map <- region_map[
    !duplicated(region_map[[marker_col]]),
    ,
    drop = FALSE
  ]
  
  markers <- region_map[[marker_col]]
  
  if (length(markers) < 2) {
    stop(
      "Fewer than two markers remain in the selected region."
    )
  }
  
  #--------------------------------------------------#
  # 5. Prepare the regional LD data
  #--------------------------------------------------#
  
  r2_region <- r2[
    markers,
    markers,
    drop = FALSE
  ]
  
  marker_positions <- stats::setNames(
    region_map[[pos_col]],
    markers
  )
  
  ld_data <- as.data.frame(
    as.table(r2_region),
    stringsAsFactors = FALSE
  )
  
  names(ld_data) <- c(
    "marker1",
    "marker2",
    "r2"
  )
  
  ld_data$marker1 <- as.character(
    ld_data$marker1
  )
  
  ld_data$marker2 <- as.character(
    ld_data$marker2
  )
  
  ld_data$position1 <- marker_positions[
    ld_data$marker1
  ]
  
  ld_data$position2 <- marker_positions[
    ld_data$marker2
  ]
  
  ld_data$index1 <- match(
    ld_data$marker1,
    markers
  )
  
  ld_data$index2 <- match(
    ld_data$marker2,
    markers
  )
  
  # Keep one triangular half and remove self-comparisons.
  ld_data <- ld_data[
    ld_data$index2 < ld_data$index1 &
      !is.na(ld_data$r2),
    ,
    drop = FALSE
  ]
  
  ld_data$marker_distance <- abs(
    ld_data$position1 -
      ld_data$position2
  )
  
  if (!is.null(max_ld_distance)) {
    
    ld_data <- ld_data[
      ld_data$marker_distance <= max_ld_distance,
      ,
      drop = FALSE
    ]
  }
  
  if (nrow(ld_data) == 0) {
    stop(
      "No pairwise LD comparisons remain after filtering."
    )
  }
  
  # Coordinates for the diamond-shaped LD plot.
  ld_data$x <- (
    ld_data$position1 +
      ld_data$position2
  ) / 2
  
  # Vertical position represents half the physical
  # distance between each pair of markers.
  ld_data$y <- ld_data$marker_distance / 2
  
  #--------------------------------------------------#
  # 6. Extract genes from the GFF table
  #--------------------------------------------------#
  
  genes <- gff_table[
    gff_table$chr_short == chromosome_short &
      as.character(gff_table$type) == "gene" &
      gff_table$start <= region_end &
      gff_table$end >= region_start,
    ,
    drop = FALSE
  ]
  
  genes <- genes[
    order(genes$start),
    ,
    drop = FALSE
  ]
  
  if (nrow(genes) == 0) {
    warning(
      "No genes were found on ",
      chromosome_short,
      " between ",
      format(region_start, big.mark = ","),
      " and ",
      format(region_end, big.mark = ","),
      " bp."
    )
  }
  
  #--------------------------------------------------#
  # 7. Prepare gene labels
  #--------------------------------------------------#
  
  if ("Name" %in% names(genes)) {
    
    genes$gene_label <- as.character(
      genes$Name
    )
    
  } else if ("ID" %in% names(genes)) {
    
    genes$gene_label <- sub(
      "^.*\\.",
      "",
      as.character(genes$ID)
    )
    
  } else {
    
    genes$gene_label <- paste0(
      "Gene_",
      seq_len(nrow(genes))
    )
  }
  
  genes$gene_midpoint <- (
    genes$start +
      genes$end
  ) / 2
  
  #--------------------------------------------------#
  # 8. Identify relevant genes
  #--------------------------------------------------#
  
  annotation_columns <- intersect(
    c(
      "Name",
      "ID",
      "Note",
      "description",
      "Description",
      "Dbxref",
      "Ontology_term"
    ),
    names(genes)
  )
  
  if (length(annotation_columns) > 0) {
    
    annotation_text <- apply(
      genes[, annotation_columns, drop = FALSE],
      1,
      function(x) {
        paste(
          x[!is.na(x)],
          collapse = " "
        )
      }
    )
    
  } else {
    
    annotation_text <- genes$gene_label
  }
  
  genes$is_gene <- grepl(
    gene_pattern,
    annotation_text,
    ignore.case = TRUE
  )
  candidate_label <- paste(gene_name, "gene")
  
  genes$gene_class <- ifelse(
    genes$is_gene,
    candidate_label,
    "Other gene"
  )
  
  #--------------------------------------------------#
  # 9. Assign nonoverlapping gene-bar tracks
  #--------------------------------------------------#
  
  assign_gene_tracks <- function(
    starts,
    ends,
    padding
  ) {
    
    tracks <- integer(length(starts))
    track_ends <- numeric(0)
    
    for (i in seq_along(starts)) {
      
      placed <- FALSE
      
      if (length(track_ends) > 0) {
        
        for (track_number in seq_along(track_ends)) {
          
          if (
            starts[i] >
            track_ends[track_number] + padding
          ) {
            
            tracks[i] <- track_number
            track_ends[track_number] <- ends[i]
            placed <- TRUE
            break
          }
        }
      }
      
      if (!placed) {
        
        tracks[i] <- length(track_ends) + 1
        
        track_ends <- c(
          track_ends,
          ends[i]
        )
      }
    }
    
    tracks
  }
  
  region_width <- region_end - region_start
  
  if (nrow(genes) > 0) {
    
    genes$track <- assign_gene_tracks(
      starts = genes$start,
      ends = genes$end,
      padding = region_width * 0.012
    )
  }
  
  max_gene_track <- if (nrow(genes) > 0) {
    max(genes$track)
  } else {
    1
  }
  
  #--------------------------------------------------#
  # 10. Assign nonoverlapping gene label lanes
  #--------------------------------------------------#
  
  assign_label_lanes <- function(
    label_midpoints,
    label_text,
    region_width,
    character_width_fraction
  ) {
    
    if (length(label_midpoints) == 0) {
      return(integer(0))
    }
    
    order_index <- order(
      label_midpoints
    )
    
    ordered_midpoints <- label_midpoints[
      order_index
    ]
    
    ordered_text <- label_text[
      order_index
    ]
    
    # Estimate the horizontal width occupied by each label.
    # Labels are horizontal, so longer gene names receive
    # wider collision intervals.
    estimated_width <- pmax(
      nchar(ordered_text),
      8
    ) *
      region_width *
      character_width_fraction
    
    label_left <- ordered_midpoints -
      estimated_width / 2
    
    label_right <- ordered_midpoints +
      estimated_width / 2
    
    lane_assignment <- integer(
      length(ordered_midpoints)
    )
    
    lane_right_edge <- numeric(0)
    
    for (i in seq_along(ordered_midpoints)) {
      
      placed <- FALSE
      
      if (length(lane_right_edge) > 0) {
        
        for (lane_number in seq_along(lane_right_edge)) {
          
          if (
            label_left[i] >
            lane_right_edge[lane_number]
          ) {
            
            lane_assignment[i] <- lane_number
            lane_right_edge[lane_number] <- label_right[i]
            placed <- TRUE
            break
          }
        }
      }
      
      if (!placed) {
        
        lane_assignment[i] <- length(
          lane_right_edge
        ) + 1
        
        lane_right_edge <- c(
          lane_right_edge,
          label_right[i]
        )
      }
    }
    
    result <- integer(
      length(label_midpoints)
    )
    
    result[order_index] <- lane_assignment
    
    result
  }
  
  genes$label_lane <- NA_integer_
  genes$label_y <- NA_real_
  
  if (
    nrow(genes) > 0 &&
    any(genes$is_gene)
  ) {
    
    gene_indices <- which(
      genes$is_gene
    )
    
    genes$label_lane[gene_indices] <-
      assign_label_lanes(
        label_midpoints = genes$gene_midpoint[
          gene_indices
        ],
        label_text = genes$gene_label[
          gene_indices
        ],
        region_width = region_width,
        character_width_fraction = label_char_width
      )
    
    genes$label_y[gene_indices] <-
      max_gene_track +
      0.60 +
      (
        genes$label_lane[gene_indices] - 1
      ) *
      label_lane_spacing
  }
  
  max_label_y <- if (
    nrow(genes) > 0 &&
    any(genes$is_gene)
  ) {
    
    max(
      genes$label_y[
        genes$is_gene
      ],
      na.rm = TRUE
    )
    
  } else {
    
    max_gene_track + 0.60
  }
  
  #--------------------------------------------------#
  # 11. Prepare significant SNP labels
  #--------------------------------------------------#
  
  significant_map <- region_map[
    FALSE,
    ,
    drop = FALSE
  ]
  
  if (!is.null(significant_snps)) {
    
    significant_snps <- as.character(
      significant_snps
    )
    
    significant_map <- region_map[
      region_map[[marker_col]] %in% significant_snps,
      ,
      drop = FALSE
    ]
    
    significant_map$snp_order <- match(
      significant_map[[marker_col]],
      significant_snps
    )
    
    significant_map <- significant_map[
      order(significant_map$snp_order),
      ,
      drop = FALSE
    ]
    
    significant_map$snp_short <- paste0(
      "S",
      significant_map$snp_order
    )
  }
  
  # Position significant SNP labels above the gene bars,
  # but below the first gene-name lane when possible.
  snp_label_y <- max_gene_track + 0.25
  
  #--------------------------------------------------#
  # 12. Create the shared physical-position scale
  #--------------------------------------------------#
  
  shared_x_scale <- ggplot2::scale_x_continuous(
    limits = c(
      region_start,
      region_end
    ),
    labels = scales::label_number(
      scale = 1e-6,
      suffix = " Mb",
      accuracy = 0.01
    ),
    expand = ggplot2::expansion(
      mult = c(0.02, 0.05)
    )
  )
  
  #--------------------------------------------------#
  # 13. Create the gene panel
  #--------------------------------------------------#
  
  gene_plot <- ggplot2::ggplot()
  
  if (nrow(genes) > 0) {
    
    # Draw all genes as horizontal bars.
    gene_plot <- gene_plot +
      ggplot2::geom_segment(
        data = genes,
        ggplot2::aes(
          x = .data[['start']],
          xend = .data[['end']],
          y = .data[['track']],
          yend = .data[['track']],
          color =.data [['gene_class']]
        ),
        linewidth = 3.7,
        lineend = "butt"
      )
    
    # Label relavant genes using fixed nonoverlapping lanes.
    # No connector lines are drawn.
    if (
      label_candidate_genes &&
      any(genes$is_gene)
    ) {
      
      gene_plot <- gene_plot +
        ggplot2::geom_text(
          data = genes[
            genes$is_gene,
            ,
            drop = FALSE
          ],
          ggplot2::aes(
            x = .data[['gene_midpoint']],
            y = .data[['label_y']],
            label = .data[['gene_label']]
          ),
          angle = 0,
          hjust = 0.5,
          vjust = 0.5,
          size = 3.0,
          fontface = "bold",
          color = "firebrick4",
          check_overlap = FALSE
        )
    }
  }
  
  # Add significant SNP lines and short labels.
  if (nrow(significant_map) > 0) {
    
    gene_plot <- gene_plot +
      ggplot2::geom_vline(
        data = significant_map,
        ggplot2::aes(
          xintercept = .data[[pos_col]]
        ),
        color = "black",
        linetype = "dashed",
        linewidth = 0.45
      ) +
      ggplot2::geom_label(
        data = significant_map,
        ggplot2::aes(
          x = .data[[pos_col]],
          y = snp_label_y,
          label = .data[['snp_short']]
        ),
        size = 3.0,
        fontface = "bold",
        fill = "white",
        linewidth = 0.2,
        label.padding = grid::unit(
          0.08,
          "lines"
        )
      )
  }
  
  gene_panel_top <- max(
    max_label_y + 0.45,
    snp_label_y + 0.45
  )
  candidate_label <- paste(gene_name, "gene")
  
  gene_plot <- gene_plot +
    ggplot2::scale_color_manual(
      values = c(
        setNames("firebrick3", candidate_label),
        "Other gene" = "grey55"
      ),
      name = "Gene class"
    )+
    shared_x_scale +
    ggplot2::scale_y_continuous(
      limits = c(
        0.5,
        gene_panel_top
      ),
      breaks = NULL,
      expand = ggplot2::expansion(
        mult = c(0.02, 0.03)
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = c(
        region_start,
        region_end
      ),
      ylim = c(
        0.5,
        gene_panel_top
      ),
      clip = "off"
    ) +
    ggplot2::labs(
      title = paste0(
        chromosome_short,
        " candidate-gene region"
      ),
      subtitle = paste(gene_name, "genes are labeled"),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12
      ),
      plot.subtitle = ggplot2::element_text(
        size = 9
      ),
      plot.margin = ggplot2::margin(
        t = 5,
        r = 65,
        b = 5,
        l = 5
      ),
      legend.position = "right"
    )
  
  #--------------------------------------------------#
  # 14. Create the LD panel
  #--------------------------------------------------#
  
  ld_plot <- ggplot2::ggplot(
    ld_data,
    ggplot2::aes(
      x = .data[['x']],
      y = .data[['y']],
      fill = .data[['r2']]
    )
  ) +
    ggplot2::geom_point(
      shape = 23,
      size = 5,
      stroke = 0.2,
      color = "grey65"
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "red3",
      limits = c(0, 1),
      na.value = "grey90",
      name = expression(r^2)
    ) +
    shared_x_scale +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(
        scale = 1e-3,
        suffix = " kb",
        accuracy = 1
      ),
      breaks = scales::breaks_pretty(
        n = 5
      ),
      expand = ggplot2::expansion(
        mult = c(0.12, 0.05)
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = c(
        region_start,
        region_end
      ),
      clip = "off"
    ) +
    ggplot2::labs(
      title = expression(
        "Pairwise linkage disequilibrium (" *
          r^2 *
          ")"
      ),
      subtitle = paste0(
        "Diagonal self-comparisons removed; ",
        length(markers),
        " markers"
      ),
      x = paste0(
        chromosome_short,
        " position"
      ),
      y = "Pairwise SNP distance / 2"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12
      ),
      plot.subtitle = ggplot2::element_text(
        size = 9
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(
          r = 8
        )
      ),
      plot.margin = ggplot2::margin(
        t = 5,
        r = 65,
        b = 8,
        l = 5
      ),
      legend.position = "right"
    )
  
  if (nrow(significant_map) > 0) {
    
    ld_plot <- ld_plot +
      ggplot2::geom_vline(
        data = significant_map,
        ggplot2::aes(
          xintercept = .data[[pos_col]]
        ),
        color = "black",
        linetype = "dashed",
        linewidth = 0.35
      )
  }
  
  #--------------------------------------------------#
  # 15. Prepare the significant-SNP key
  #--------------------------------------------------#
  
  if (nrow(significant_map) > 0) {
    
    snp_key_text <- paste(
      paste0(
        significant_map$snp_short,
        " = ",
        significant_map[[marker_col]]
      ),
      collapse = "; "
    )
    
  } else {
    
    snp_key_text <- NULL
  }
  
  #--------------------------------------------------#
  # 16. Combine the panels
  #--------------------------------------------------#
  
  number_label_lanes <- if (
    nrow(genes) > 0 &&
    any(genes$is_gene)
  ) {
    
    max(
      genes$label_lane[
        genes$is_gene
      ],
      na.rm = TRUE
    )
    
  } else {
    
    1
  }
  
  gene_panel_height <- max(
    1.45,
    1.15 + number_label_lanes * 0.20
  )
  
  combined_plot <- (
    gene_plot /
      ld_plot
  ) +
    patchwork::plot_layout(
      heights = c(
        gene_panel_height,
        2.4
      ),
      guides = "collect"
    ) +
    patchwork::plot_annotation(
      tag_levels = "a",
      caption = snp_key_text,
      theme = ggplot2::theme(
        plot.tag = ggplot2::element_text(
          face = "bold",
          size = 13
        ),
        plot.caption = ggplot2::element_text(
          hjust = 0,
          size = 9
        )
      )
    )
  
  print(combined_plot)
  
  #--------------------------------------------------#
  # 17. Return the plot and supporting data
  #--------------------------------------------------#
  
  invisible(
    list(
      plot = combined_plot,
      gene_plot = gene_plot,
      ld_plot = ld_plot,
      genes = genes,
      candidate_genes = genes[
        genes$is_gene,
        ,
        drop = FALSE
      ],
      marker_map = region_map,
      significant_map = significant_map,
      r2 = r2_region,
      ld_data = ld_data,
      chromosome = chromosome_short,
      region_start = region_start,
      region_end = region_end
    )
  )
}
