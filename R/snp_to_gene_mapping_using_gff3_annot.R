#' Map SNPs to Nearby Genes Using GFF3 Genome Annotations
#'
#' Maps a set of target SNPs to nearby gene models using genomic coordinates
#' from a GFF3 annotation file. For each SNP, the function searches upstream
#' and downstream within a user-defined physical window and identifies all
#' genes whose genomic intervals overlap that window.
#'
#' The function is designed for post-GWAS candidate-gene analysis. It first
#' synchronizes SNP identifiers and genomic coordinates with the supplied
#' physical map, imports the GFF3 annotation using \pkg{rtracklayer}, and
#' extracts only features classified as genes. Gene coordinates and selected
#' annotation metadata are then used to identify genes located within the
#' specified distance of each target SNP.
#'
#' Gene position relative to each SNP is classified as \code{"Upstream"},
#' \code{"Downstream"}, or \code{"Overlapping"} based on the position of the
#' gene midpoint relative to the SNP coordinate.
#'
#' This function currently expects chromosome identifiers that contain the
#' soybean chromosome naming convention \code{"Gm01"}--\code{"Gm20"}. It is
#' therefore primarily intended for soybean GFF3 annotations unless the
#' chromosome-extraction step is modified for another genome.
#'
#' @param snp_list Character vector of target SNP identifiers to query.
#'   Typically, this is a vector containing significant SNPs identified from
#'   a GWAS or other marker-trait association analysis.
#' @param gff_path Character string specifying the path to a local GFF3
#'   annotation file. The file should contain genomic coordinates and gene
#'   annotations. For the soybean Williams 82 reference genome, annotations
#'   can be obtained from the SoyBase genome annotation repository:
#'   \url{https://data.soybase.org/Glycine/max/annotations/Wm82.gnm1.ann1.DvBy/}.
#' @param map_data A data frame containing physical marker coordinates.
#'   It must contain columns named \code{SNP}, \code{Chromosome}, and
#'   \code{Position}. \code{SNP} must contain marker identifiers matching
#'   those supplied in \code{snp_list}.
#' @param window_bp Numeric value specifying the upstream and downstream
#'   physical search distance around each SNP, in base pairs. A gene is
#'   retained when its genomic interval overlaps the SNP-centered window.
#'   Default is \code{50000} (50 kb).
#'
#' @return A tibble containing one row for each SNP-gene association found
#'   within the specified physical window. The returned table contains:
#'   \itemize{
#'     \item \code{query_snp}: Target SNP identifier.
#'     \item \code{snp_chr}: Chromosome containing the SNP.
#'     \item \code{snp_pos}: Physical position of the SNP in base pairs.
#'     \item \code{gene_id}: Gene identifier extracted from the GFF3 annotation.
#'     \item \code{gene_symbol}: Gene name or symbol extracted from the GFF3
#'       \code{Name} attribute, when available.
#'     \item \code{seqnames}: Original chromosome/sequence identifier from
#'       the GFF3 file.
#'     \item \code{start}: Gene start coordinate.
#'     \item \code{end}: Gene end coordinate.
#'     \item \code{strand}: Gene strand orientation.
#'     \item \code{type}: GFF3 feature type; retained as \code{"gene"}.
#'     \item \code{distance_to_snp}: Distance in base pairs between the SNP
#'       position and the midpoint of the gene. Positive values indicate that
#'       the gene midpoint is downstream of the SNP, while negative values
#'       indicate that it is upstream.
#'     \item \code{orientation}: Relative position of the gene with respect
#'       to the SNP: \code{"Upstream"}, \code{"Downstream"}, or
#'       \code{"Overlapping"}.
#'     \item \code{description}: Functional or descriptive annotation
#'       extracted from the GFF3 \code{Note} field, when available.
#'   }
#'
#' @details
#' The mapping procedure consists of four main steps:
#'
#' \enumerate{
#'   \item Target SNP identifiers are cleaned and matched against the
#'   supplied physical map.
#'
#'   \item The GFF3 file is imported and filtered to retain only features
#'   with \code{type == "gene"}. Selected attributes, including gene ID,
#'   gene name, and description, are extracted into a lightweight data frame.
#'
#'   \item SNPs and genes are matched by chromosome, followed by an interval
#'   overlap test. A gene is retained when its start coordinate is before
#'   the end of the SNP-centered search window and its end coordinate is
#'   after the beginning of that window.
#'
#'   \item The distance between each SNP and the midpoint of each nearby
#'   gene is calculated and used to classify the gene as upstream,
#'   downstream, or overlapping.
#' }
#'
#' The function does not modify the original genotype or map data. It returns
#' only SNP-gene relationships detected within the specified physical window.
#'
#' @section Interpretation:
#' A gene is considered "Overlapping" when its midpoint is exactly equal to
#' the SNP position. This classification is based on gene midpoint rather
#' than gene interval overlap. A gene can therefore physically overlap the
#' SNP while still being classified as upstream or downstream if its midpoint
#' lies on one side of the SNP.
#'
#' @section Genome compatibility:
#' Chromosome identifiers are standardized by extracting patterns matching
#' \code{Gm\\d+}. Consequently, this implementation is optimized for soybean
#' chromosome naming conventions such as \code{Gm01}, \code{Gm02}, and
#' \code{Gm20}. For other crops, the chromosome-standardization step should
#' be adapted to the naming convention used by the corresponding GFF3 file.
#'
#' @export
#'
#' @importFrom rtracklayer import
#' @importFrom dplyr mutate filter select arrange tibble case_when inner_join
#' @importFrom stringr str_extract
#'
#' @examples
#' \dontrun{
#' library(LDdecay)
#'
#' # Significant SNPs identified from a GWAS
#' test_snps <- c(
#'   "ss715620779",
#'   "ss715620778",
#'   "ss715620770"
#' )
#'
#' # Map SNPs to genes within 50 kb
#' test_mapping <- snp_to_gene_mapping_using_gff3_annot(
#'   snp_list = test_snps,
#'   gff_path = "data/glyma.Wm82.gnm1.ann1.DvBy.gene_models_main.gff3",
#'   map_data = map,
#'   window_bp = 50000
#' )
#'
#' # Inspect the SNP-gene associations
#' head(test_mapping)
#'
#' # View genes surrounding each target SNP
#' print(test_mapping)
#' }
#'
snp_to_gene_mapping_using_gff3_annot <- function(
    snp_list,
    gff_path,
    map_data,
    window_bp = 50000
) {
  
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
  
  if (nrow(gff_genes) == 0) {
    stop("Error: No structural feature tags labeled 'gene' found in the GFF3 file.")
  }
  
  message("-> Stripping and converting structural elements to minimal tables...")
  
  # Safe extraction of attributes without calling un-vectorized df inflation
  gff_meta <- S4Vectors::ncol(gff_genes)
  
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
