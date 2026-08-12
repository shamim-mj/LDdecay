# Plot Linkage Disequilibrium Heatmap Around a Focal Variant

Calculates pairwise linkage disequilibrium (LD) as squared Pearson
correlations (\\r^2\\) among variants surrounding a focal SNP and
visualizes the resulting LD matrix as a heatmap. Variants are restricted
to the same chromosome as the focal SNP, filtered by a user-defined
physical window, and screened using a minimum minor allele frequency
(MAF) threshold.

## Usage

``` r
plot_single_snp_ld(
  geno,
  map,
  focal_snp,
  window_bp = 1e+05,
  marker_col = "SNP",
  chr_col = "Chromosome",
  pos_col = "Position",
  min_maf = 0.05
)
```

## Arguments

- geno:

  A data frame or matrix containing genotype dosages, with rows
  representing samples and columns representing variants. If present, a
  column named `"taxa"` is automatically removed before LD calculations.
  Genotype values should be numerically coded, such as 0, 1, and 2.

- map:

  A data frame containing physical marker information. It must contain
  columns identifying the variant, chromosome, and physical position.
  The corresponding column names are specified using `marker_col`,
  `chr_col`, and `pos_col`.

- focal_snp:

  Character string specifying the focal variant around which the LD
  region is constructed. The focal SNP must be present in both the
  marker map and genotype data.

- window_bp:

  Numeric. Physical distance in base pairs upstream and downstream of
  the focal SNP to include in the analysis. Default is `100000` (100 kb
  on each side of the focal SNP).

- marker_col:

  Character string specifying the column in `map` containing variant
  identifiers. Default is `"SNP"`.

- chr_col:

  Character string specifying the column in `map` containing chromosome
  identifiers. Default is `"Chromosome"`.

- pos_col:

  Character string specifying the column in `map` containing physical
  marker positions in base pairs. Default is `"Position"`.

- min_maf:

  Numeric. Minimum minor allele frequency required for a variant to be
  retained in the LD analysis. Variants with MAF below this threshold
  are removed. Default is `0.05`.

## Value

A list containing:

- plot:

  A `ggplot` object showing the pairwise local LD matrix as a heatmap.

- r2:

  A square matrix containing pairwise squared Pearson correlation
  coefficients (\\r^2\\) among the retained variants.

- map:

  A data frame containing the marker metadata for variants retained
  after regional and MAF filtering.

- genotypes:

  A numeric matrix containing genotype dosages for the retained
  variants.

## Details

The focal SNP is highlighted with an asterisk in the heatmap axis
labels. The resulting plot provides a visual representation of local LD
structure around a significant or otherwise prioritized variant.

The function first identifies the chromosome and physical position of
the focal SNP. It then selects markers located on the same chromosome
and within `window_bp` base pairs upstream or downstream of the focal
variant.

Variants not present in the genotype matrix are removed. Minor allele
frequency is calculated from genotype dosage values as:

\$\$ AF = \frac{\mathrm{mean}(G)}{2} \$\$

and minor allele frequency is calculated as:

\$\$ MAF = \min(AF, 1 - AF). \$\$

Variants with MAF below `min_maf` are excluded before calculating
pairwise LD.

Pairwise LD is calculated as the squared Pearson correlation between
genotype dosage vectors:

\$\$ r^2 = cor(G_i, G_j)^2. \$\$

Correlations are calculated using `use = "pairwise.complete.obs"`,
allowing each pairwise comparison to use available non-missing genotype
observations.

The heatmap displays the lower triangular portion of the pairwise
correlation matrix. The focal SNP is identified by an asterisk in the
axis labels.

## Examples

``` r
if (FALSE) { # \dontrun{
library(LDdecay)

# Plot LD surrounding a focal SNP using a 100-kb window
ld_result <- plot_single_snp_ld(
  geno = geno,
  map = map,
  focal_snp = "ss715620779",
  window_bp = 100000
)

# Display the heatmap
ld_result$plot

# Inspect the pairwise LD matrix
head(ld_result$r2)

# Inspect the markers retained in the analysis
head(ld_result$map)
} # }
```
