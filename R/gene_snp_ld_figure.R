#' Plot Linkage Disequilibrium Matrix and Gene Annotations in an Area of Interest (AOI)
#'
#' Renders a publication-quality composite plot aligning a pairwise Linkage Disequilibrium (LD) 
#' correlation heatmap with structural gene models extracted from a GFF database. Highlights
#' candidate gene sets matching key functional patterns and tracks user-defined variant subsets.
#'
#' @param map A data frame containing physical variant metadata map records. Must include columns matching `SNP`, `Chr`, and `Position`.
#' @param gff_table A processed data frame containing parsed genomic feature fields (e.g., coordinates, IDs, Names, and Notes). Use the 'get_gff' function before running this code
#' @param chromosome Character string or numeric value specifying the single target chromosome to visualize (e.g., "Gm03").
#' @param focal_snp Character string identifying the central target variant identifier (e.g., `"ss715620779"`).
#' @param ld_window Numeric. The search window radius in base-pairs flanking the variant. Default is set to `360000` (360kb flanking radius).
#' @param region_start Numeric. The minimum base-pair boundary coordinate defining the active genomic window. Default is `NULL` (automatic min coordinate).
#' @param region_end Numeric. The maximum base-pair boundary coordinate defining the active genomic window. Default is `NULL` (automatic max coordinate).
#' @param marker_col Character. Name of the column containing unique variant identifiers in the `map` file. Default is `"marker"`.
#' @param chr_col Character. Name of the column containing chromosome assignments in the `map` file. Default is `"chr"`.
#' @param pos_col Character. Name of the column containing base-pair coordinates in the `map` file. Default is `"position"`.
#' @param significant_snps A character vector of prioritized or highly significant variant identifiers to visually highlight on the plot framework. Default is `NULL`.
#' @param gene_pattern Character string containing a regular expression (regex) to isolate and highlight priority target pathways in descriptions. Default targets specific ABC/PDR plant transporter annotations.
#' @param gene_name Character. Title label or text indicator assigned to the structural gene tracks matching `cyp_pattern` in plot legends. Default is `"ABC/PDR Transport"`.
#' @param label_cyp_genes Logical. If `TRUE`, includes text annotation callouts on top of structural gene models matching your target functional pattern. Default is `TRUE`.
#' @param max_ld_distance Numeric. Maximum physical distance cutoff limit (in base pairs) to bound or clip pairwise LD calculations. Default is `NULL`.
#' @param label_char_width Numeric. Text scaling modifier establishing character width sizing constraints during layout evaluations. Default is `0.010`.
#' @param label_lane_spacing Numeric. Vertical adjustment increment determining track separation margins for overlapping gene model callouts. Default is `0.55`.
#'
#' @return A composite plot object layout aligning physical coordinates, structural annotations, and decay profiles.
#' @export
#'
#' @import ggplot2
#' @importFrom dplyr filter mutate select
#' @importFrom stringr str_detect
plot_aoi_ld_region <- function(
    map,
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
    label_cyp_genes = TRUE, max_ld_distance = NULL, 
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
  
  ld_result <- plot_snp_ld(
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
  
  marker_positions <- setNames(
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
  # 8. Identify CYP450 genes
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
  
  genes$is_CYP450 <- grepl(
    gene_pattern,
    annotation_text,
    ignore.case = TRUE
  )
  
  genes$gene_class <- ifelse(
    genes$is_CYP450,
    paste(gene_name, "gene"),
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
  # 10. Assign nonoverlapping CYP450 label lanes
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
    any(genes$is_CYP450)
  ) {
    
    cyp_indices <- which(
      genes$is_CYP450
    )
    
    genes$label_lane[cyp_indices] <-
      assign_label_lanes(
        label_midpoints = genes$gene_midpoint[
          cyp_indices
        ],
        label_text = genes$gene_label[
          cyp_indices
        ],
        region_width = region_width,
        character_width_fraction = label_char_width
      )
    
    genes$label_y[cyp_indices] <-
      max_gene_track +
      0.60 +
      (
        genes$label_lane[cyp_indices] - 1
      ) *
      label_lane_spacing
  }
  
  max_label_y <- if (
    nrow(genes) > 0 &&
    any(genes$is_CYP450)
  ) {
    
    max(
      genes$label_y[
        genes$is_CYP450
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
          x = start,
          xend = end,
          y = track,
          yend = track,
          color = gene_class
        ),
        linewidth = 3.7,
        lineend = "butt"
      )
    
    # Label CYP450 genes using fixed nonoverlapping lanes.
    # No connector lines are drawn.
    if (
      label_cyp_genes &&
      any(genes$is_CYP450)
    ) {
      
      gene_plot <- gene_plot +
        ggplot2::geom_text(
          data = genes[
            genes$is_CYP450,
            ,
            drop = FALSE
          ],
          ggplot2::aes(
            x = gene_midpoint,
            y = label_y,
            label = gene_label
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
          label = snp_short
        ),
        size = 3.0,
        fontface = "bold",
        fill = "white",
        label.size = 0.2,
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
        candidate_label = "firebrick3",
        "Other gene" = "grey55"
      ),
      name = "Gene class"
    ) +
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
      x = x,
      y = y,
      fill = r2
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
    any(genes$is_CYP450)
  ) {
    
    max(
      genes$label_lane[
        genes$is_CYP450
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
      cyp_genes = genes[
        genes$is_CYP450,
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
