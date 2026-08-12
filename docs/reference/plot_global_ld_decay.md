# Calculate and Plot Global Linkage Disequilibrium Decay

Calculates genome-wide pairwise linkage disequilibrium (LD) as squared
Pearson correlations (\\r^2\\) between genotype markers within a
user-defined physical distance. The function synchronizes genotype and
marker-map data, removes unanchored scaffold markers, calculates
pairwise \\r^2\\ values within chromosomes, and generates a scatterplot
of LD against physical distance. A nonlinear least-squares (NLS)
regression curve is fitted using the model \\r^2 = 1/(1 + C d)\\, where
\\d\\ is physical distance and \\C\\ is the fitted decay parameter.

## Usage

``` r
plot_global_ld_decay(
  geno,
  map,
  max_distance_bp = 5e+05,
  max_plot_points = 1e+05,
  r2_threshold = 0.2,
  output_image_path = "Figures/ld_short_range.png"
)
```

## Arguments

- geno:

  A data frame or matrix containing genotype dosages. A column named
  `"taxa"` is required to identify samples, while the remaining columns
  should contain variant identifiers matching the `SNP` column in `map`.
  Genotype values should be coded numerically (e.g., 0, 1, and 2).

- map:

  A data frame containing physical marker information. It must contain
  the columns `"SNP"`, `"Chromosome"`, and `"Position"`, corresponding
  to variant identifier, chromosome, and physical position in base
  pairs, respectively.

- max_distance_bp:

  Numeric. Maximum physical distance in base pairs between two markers
  for which pairwise LD is calculated. Default is `500000` (500 kb).

- max_plot_points:

  Numeric. Maximum number of pairwise LD observations randomly sampled
  for visualization. This limits the number of points rendered in the
  plot and helps reduce memory and plotting requirements. Default is
  `100000`.

- r2_threshold:

  Numeric. Value used to draw a horizontal reference line on the plot,
  representing a user-defined LD threshold. Default is `0.2`.

- output_image_path:

  Character string specifying the path, including file extension, where
  the generated plot should be saved. If `NULL`, the plot is not saved
  to a file. Default is `"Figures/ld_short_range.png"`.

## Value

A list containing:

- plot:

  A `ggplot` object showing pairwise \\r^2\\ against physical distance,
  with an NLS LD-decay curve and user-defined \\r^2\\ reference line.

- ld_data:

  A data frame containing the calculated pairwise LD observations after
  removal of missing values and self-comparisons. Columns include `d`,
  the physical distance between marker pairs, and `r2`, the
  corresponding squared correlation coefficient.

## Details

To reduce memory and rendering requirements, a maximum number of
pairwise observations can be randomly sampled for visualization. The
returned `ld_data` object contains the complete set of calculated
pairwise observations retained after filtering, whereas the plot is
generated from the sampled observations.

Markers classified as unanchored scaffolds are removed before LD
calculations. Only variants present in both the genotype and map
datasets are retained. Pairwise LD is calculated only between markers
located on the same chromosome and within `max_distance_bp` base pairs.

The fitted decay model is:

\$\$ r^2 = \frac{1}{1 + C d} \$\$

where \\d\\ is physical distance and \\C\\ is estimated using nonlinear
least squares. The regression curve is fitted to a random sample of at
most `max_plot_points` pairwise observations. The complete filtered
pairwise dataset is returned in `ld_data`.

## Examples

``` r
if (FALSE) { # \dontrun{
library(LDdecay)

# Calculate genome-wide LD decay
ld_result <- plot_global_ld_decay(
  geno = geno,
  map = map,
  max_distance_bp = 500000,
  max_plot_points = 100000,
  r2_threshold = 0.2
)

# Display the LD-decay plot
ld_result$plot

# Inspect the calculated pairwise LD data
head(ld_result$ld_data)

# Save the plot to a custom location
ld_result <- plot_global_ld_decay(
  geno = geno,
  map = map,
  output_image_path = "Figures/global_ld_decay.png"
)
} # }
```
