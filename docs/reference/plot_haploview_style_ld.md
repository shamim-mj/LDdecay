# Generate a Haploview-Style Inverted Triangle LD Heatmap

Calculates pairwise linkage disequilibrium (LD) as squared Pearson
correlations (\\r^2\\) among a selected set of variants and renders the
resulting matrix as a Haploview-style LD heatmap. The heatmap can be
displayed in an inverted-triangle configuration and is aligned with the
physical positions of the selected variants.

## Usage

``` r
plot_haploview_style_ld(
  geno,
  map,
  snp_subset = NULL,
  title = "Pairwise LD",
  color_palette = grDevices::heat.colors(20),
  flip = TRUE
)
```

## Arguments

- geno:

  A data frame or matrix containing genotype dosages. Columns
  corresponding to the variants specified in `snp_subset` must be
  present in the genotype data. Genotypes should be coded numerically
  (e.g., 0, 1, and 2).

- map:

  A data frame containing physical marker information. It must contain
  the columns `"SNP"` and `"Position"`, corresponding to variant
  identifiers and physical positions in base pairs, respectively. A
  `"Chromosome"` column may also be present for compatibility with
  standard marker-map files.

- snp_subset:

  Specifies the variants to include in the heatmap. The argument can be:

  - a character vector of SNP identifiers;

  - a numeric vector of row indices in `map`; or

  - a single character range such as `"1:60"`.

  If `NULL`, the first 60 variants in `map` are selected.

- title:

  Character string specifying the title displayed on the heatmap.
  Default is `"Pairwise LD"`.

- color_palette:

  Character vector defining the colors used to represent LD values. The
  default is `grDevices::heat.colors(20)`.

- flip:

  Logical. If `TRUE`, the LDheatmap is displayed in the
  inverted-triangle orientation. Default is `TRUE`.

## Value

An object returned by
[`LDheatmap::LDheatmap`](https://github.com/mcneney/LDheatmap/reference/LDheatmap.html),
containing the rendered LD heatmap and associated grid graphics
components.

## Details

The function supports several ways of specifying the variants to
display, including explicit SNP identifiers, numeric row indices, and
character ranges such as `"1:60"`. If no subset is provided, the first
60 variants in the map are selected automatically.

The function first identifies the requested variants from the marker map
and synchronizes them with the genotype matrix. Pairwise LD is
calculated as:

\$\$ r^2 = cor(G_i, G_j)^2 \$\$

where \\G_i\\ and \\G_j\\ are genotype dosage vectors for two variants.
Pairwise correlations are calculated using
`use = "pairwise.complete.obs"`, allowing available non-missing genotype
observations to be used for each pair.

The resulting square \\r^2\\ matrix is supplied to
[`LDheatmap::LDheatmap`](https://github.com/mcneney/LDheatmap/reference/LDheatmap.html),
together with the physical positions and variant identifiers. The
function therefore provides a convenient way to visualize local LD
structure for a selected group of variants.

## Examples

``` r
if (FALSE) { # \dontrun{
library(LDdecay)

# Plot the first 60 variants in the marker map
ld_plot <- plot_haploview_style_ld(
  geno = geno,
  map = map
)

# Plot variants 1 through 30
ld_plot <- plot_haploview_style_ld(
  geno = geno,
  map = map,
  snp_subset = "1:30"
)

# Plot a specific set of SNPs
selected_snps <- c(
  "ss715620779",
  "ss715620801",
  "ss715620845"
)

ld_plot <- plot_haploview_style_ld(
  geno = geno,
  map = map,
  snp_subset = selected_snps,
  title = "Candidate SNP LD"
)

# Display the heatmap
ld_plot
} # }
```
