# Calculate and Plot Binned Empirical Linkage Disequilibrium Decay

Calculates pairwise linkage disequilibrium (LD) as squared Pearson
correlation (\\r^2\\) between genotype markers and summarizes LD across
physical-distance intervals. The resulting empirical LD decay profile is
visualized as the mean \\r^2\\ across user-defined distance bins.

## Usage

``` r
plot_binned_ld_global_decay(
  geno,
  map,
  bin_size_bp = 10000,
  max_distance_bp = 5e+05,
  r2_threshold = 0.2,
  output_image_path = "Figures/ld_binned_decay.png"
)
```

## Arguments

- geno:

  A data frame or matrix containing genotype dosages. The data must
  contain a column named `"taxa"` identifying individual samples, with
  all remaining columns corresponding to variant identifiers. Genotype
  values should be numeric dosage values such as 0, 1, and 2.

- map:

  A data frame containing physical coordinates for the variants. The
  data frame must contain columns named `"SNP"`, `"Chromosome"`, and
  `"Position"`, corresponding to variant identifiers, chromosome
  assignments, and physical positions in base pairs, respectively.

- bin_size_bp:

  Numeric value specifying the width of the physical-distance bins used
  to summarize pairwise LD. The value is expressed in base pairs.
  Default is `10000` (10 kb).

- max_distance_bp:

  Numeric value specifying the maximum physical distance over which
  pairwise LD is calculated and displayed. The value is expressed in
  base pairs. Default is `500000` (500 kb).

- r2_threshold:

  Numeric value specifying the \\r^2\\ threshold displayed as a
  horizontal dashed reference line on the plot. Default is `0.2`.

- output_image_path:

  Character string specifying the file path at which the generated plot
  should be saved. The file extension determines the output format
  supported by
  [`ggsave`](https://ggplot2.tidyverse.org/reference/ggsave.html). Set
  to `NULL` to return the plot without saving an image. Default is
  `"Figures/ld_binned_decay.png"`.

## Value

A list containing two elements:

- plot:

  A `ggplot` object containing the empirical LD decay profile.

- binned_data:

  A data frame containing the summarized LD values for each
  physical-distance bin. Columns include `bin`, the lower boundary of
  the distance bin in base pairs; `mean_r2`, the mean pairwise \\r^2\\;
  `median_r2`, the median pairwise \\r^2\\; and `count`, the number of
  pairwise comparisons contributing to the bin.

## Details

Genotype and physical-map data are first synchronized using the
intersection of marker identifiers present in both datasets. Unanchored
scaffold chromosomes are excluded, and markers are ordered by chromosome
and physical position. Pairwise \\r^2\\ values are then calculated
between markers on the same chromosome within the specified maximum
physical distance.

The function calculates LD separately within each chromosome and does
not calculate LD between markers located on different chromosomes.
Markers located on chromosomes whose names begin with `"scaffold"` are
excluded before analysis.

Only markers present in both `geno` and `map` are retained. Markers are
ordered by chromosome and physical position before pairwise comparisons
are performed.

For each marker, the function identifies downstream markers on the same
chromosome within `max_distance_bp`. Pairwise LD is calculated as the
squared Pearson correlation between genotype dosage vectors using
pairwise-complete observations.

Pairwise LD values are grouped into physical-distance bins of width
`bin_size_bp`. For each bin, the function calculates the mean \\r^2\\,
median \\r^2\\, and number of contributing marker pairs. The resulting
mean \\r^2\\ values are used to construct the empirical LD decay curve.

A small calculation buffer is added internally to the requested maximum
distance to reduce edge effects when assigning pairwise observations to
the final distance bins.

If `output_image_path` is not `NULL`, the output directory is created
automatically when necessary and the plot is saved at 300 dpi.

## Examples

``` r
if (FALSE) { # \dontrun{
library(LDdecay)

# Calculate empirical LD decay using 10-kb distance bins
ld_decay <- plot_binned_ld_global_decay(
  geno = geno,
  map = map,
  bin_size_bp = 10000,
  max_distance_bp = 500000
)

# Display the LD decay plot
ld_decay$plot

# Inspect the summarized LD data
head(ld_decay$binned_data)

# Run the analysis without saving an image
ld_decay <- plot_binned_ld_global_decay(
  geno = geno,
  map = map,
  output_image_path = NULL
)
} # }
```
