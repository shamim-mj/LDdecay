# Plot Linkage Disequilibrium and Gene Annotations for a Genomic Region

Generates a publication-quality composite visualization of linkage
disequilibrium (LD) and gene annotations surrounding a focal variant.
The function calculates SNP-specific pairwise LD around the focal SNP,
extracts gene models from a GFF3-derived annotation table, identifies
candidate genes using a user-defined regular expression, and combines
the results into a two-panel genomic region plot.

## Usage

``` r
plot_snp_ld_region_and_genes(
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
  label_candidate_genes = TRUE,
  max_ld_distance = NULL,
  label_char_width = 0.01,
  label_lane_spacing = 0.55
)
```

## Arguments

- map:

  A data frame containing physical marker information. It must contain
  columns identifying the marker, chromosome, and physical position. The
  corresponding column names are specified by `marker_col`, `chr_col`,
  and `pos_col`.

- geno:

  A data frame or matrix containing genotype dosages. Rows represent
  individual samples and columns represent genetic markers. A `"taxa"`
  column identifying individuals may be included and is automatically
  removed before LD calculations. Genotype values should be numeric
  dosages, typically coded as 0, 1, and 2.

- gff_table:

  A data frame containing genomic annotations parsed from a GFF3 file.
  The table must contain at least the columns `seqnames`, `start`,
  `end`, `strand`, and `type`. Gene annotations should be available as
  records with `type == "gene"`.

- chromosome:

  Character string or numeric value identifying the chromosome to
  analyze. Chromosome identifiers are standardized to soybean-style
  names such as `"Gm03"`.

- focal_snp:

  Character string identifying the focal or significant SNP around which
  SNP-specific LD is calculated. The SNP must be present in both `map`
  and `geno`.

- ld_window:

  Numeric. Physical distance in base pairs upstream and downstream of
  `focal_snp` used when calculating SNP-specific LD. Default is `360000`
  (360 kb).

- region_start:

  Numeric. Lower physical coordinate defining the genomic region
  displayed in the final plot. If `NULL`, the minimum position among
  available markers in the selected chromosome is used.

- region_end:

  Numeric. Upper physical coordinate defining the genomic region
  displayed in the final plot. If `NULL`, the maximum position among
  available markers in the selected chromosome is used.

- marker_col:

  Character string specifying the column in `map` containing unique
  marker identifiers. Default is `"SNP"`.

- chr_col:

  Character string specifying the column in `map` containing chromosome
  identifiers. Default is `"Chromosome"`.

- pos_col:

  Character string specifying the column in `map` containing physical
  base-pair positions. Default is `"Position"`.

- significant_snps:

  Character vector containing SNP identifiers that should be highlighted
  in the genomic region plot. Significant SNPs are displayed as vertical
  dashed lines and assigned short labels (e.g., `S1`, `S2`). Default is
  `NULL`.

- gene_pattern:

  Character string containing a regular expression used to identify
  candidate genes from their available annotation fields, including gene
  names, descriptions, notes, and ontology information. Matching is
  case-insensitive. The default pattern searches for ABC/PDR
  transporter-related annotations.

- gene_name:

  Character string defining the functional category name used to label
  genes matching `gene_pattern` in the plot legend. Default is
  `"ABC/PDR Transport"`.

- label_candidate_genes:

  Logical. If `TRUE`, labels genes matching `gene_pattern` above their
  corresponding gene models. Despite the historical argument name, this
  option is not restricted to cytochrome P450 genes and can be used with
  any gene annotation pattern. Default is `TRUE`.

- max_ld_distance:

  Numeric. Optional maximum physical distance in base pairs between
  marker pairs retained for the LD panel. If `NULL`, all pairwise
  comparisons within the selected region are retained. Default is
  `NULL`.

- label_char_width:

  Numeric. Scaling factor used to estimate the horizontal space occupied
  by gene labels when assigning labels to non-overlapping lanes. Larger
  values allocate more horizontal space to each label. Default is
  `0.010`.

- label_lane_spacing:

  Numeric. Vertical spacing between successive gene-label lanes used to
  reduce label overlap. Default is `0.55`.

## Value

A named list containing the composite plot and supporting genomic data:

- `plot`:

  Complete composite `patchwork` plot containing the gene annotation and
  LD panels.

- `gene_plot`:

  The gene annotation panel as a `ggplot` object.

- `ld_plot`:

  The LD triangular panel as a `ggplot` object.

- `genes`:

  All gene models overlapping the selected genomic region, including
  gene classes, tracks, and label positions.

- `candidate_genes`:

  Subset of `genes` matching `gene_pattern`. The element name is
  retained for backward compatibility and is not restricted to CYP
  genes.

- `marker_map`:

  Marker map records retained within the selected genomic region.

- `significant_map`:

  Map records corresponding to `significant_snps`.

- `r2`:

  Square pairwise LD correlation matrix for markers retained in the
  selected region.

- `ld_data`:

  Long-format pairwise LD data containing marker identities, physical
  positions, pairwise \\r^2\\, marker distance, and plotting
  coordinates.

- `chromosome`:

  Standardized chromosome identifier used for the analysis.

- `region_start`:

  Lower physical coordinate of the plotted region.

- `region_end`:

  Upper physical coordinate of the plotted region.

## Details

The upper panel displays gene models positioned according to their
physical genomic coordinates. Genes matching `gene_pattern` are
highlighted and can be labeled automatically. Significant SNPs can also
be displayed as vertical reference lines with short labels. The lower
panel displays pairwise LD as an inverted triangular representation,
where the x-axis corresponds to physical genomic position and the y-axis
represents half the physical distance between marker pairs.

This function is designed for investigating candidate genomic regions
surrounding significant GWAS markers and can be used to evaluate whether
significant variants are located within or near genes with potentially
relevant functional annotations.

The function performs the following workflow:

1.  Calculates SNP-specific pairwise LD around the focal SNP using
    [`plot_single_snp_ld`](plot_single_snp_ld.md).

2.  Validates the marker map and GFF3-derived annotation table.

3.  Standardizes chromosome identifiers to the `Gm01`, `Gm02`, ...,
    `Gm20` format.

4.  Selects markers located on the requested chromosome and within the
    specified genomic region.

5.  Converts the pairwise LD matrix into long format and retains one
    triangular half of the pairwise comparisons.

6.  Extracts gene models overlapping the selected genomic region.

7.  Searches gene annotation fields using `gene_pattern` to identify
    candidate genes.

8.  Assigns genes to non-overlapping horizontal tracks.

9.  Assigns candidate-gene labels to non-overlapping vertical lanes.

10. Adds user-defined significant SNPs as reference lines and labels.

11. Constructs a shared physical-position scale between the gene and LD
    panels.

12. Combines the gene and LD panels using `patchwork`.

Candidate genes are identified using a regular expression applied to
available annotation fields such as `Name`, `ID`, `Note`, `Description`,
`Dbxref`, and `Ontology_term`. Therefore, users can adapt the function
to investigate different gene families or biological pathways by
changing `gene_pattern` and `gene_name`.

The LD panel uses the physical midpoint of each marker pair as its
horizontal coordinate and half of the physical marker distance as its
vertical coordinate. This produces a triangular representation similar
to commonly used LD-region visualizations.

## See also

[`plot_single_snp_ld`](plot_single_snp_ld.md), [`get_gff`](get_gff.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(LDdecay)

# Define the focal genomic region
focal_snp <- "ss715620779"
chromosome <- "Gm03"

# Define candidate-gene annotation pattern
gene_pattern <- "ABC-2/plant PDR ABC transporter|ABC|PDR transporter"
gene_name <- "ABC/PDR Transport"

# Generate LD and gene annotation plot
result <- plot_snp_ld_region_and_genes(
  map = map,
  geno = geno,
  gff_table = gff_table,
  chromosome = chromosome,
  focal_snp = focal_snp,
  ld_window = 360000,
  gene_pattern = gene_pattern,
  gene_name = gene_name
)

# Display the composite plot
result$plot

# View candidate genes
result$candidate_genes

# View pairwise LD data
head(result$ld_data)
} # }
```
