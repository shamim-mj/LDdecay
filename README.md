# LDdecay

Advanced tools for linkage disequilibrium decay analysis, regional LD visualization, and SNP-to-gene mapping.

[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue)](https://shamim-mj.github.io/LDdecay/)
[![GitHub last commit](https://img.shields.io/github/last-commit/shamim-mj/LDdecay)](https://github.com/shamim-mj/LDdecay)
[![License](https://img.shields.io/github/license/shamim-mj/LDdecay)](https://github.com/shamim-mj/LDdecay)

`LDdecay` is an R package for linkage disequilibrium (LD) analysis and visualization. The package provides tools for calculating genome-wide LD decay, generating empirical binned LD-decay curves, visualizing regional LD patterns, mapping SNPs to nearby genomic features, and integrating LD information with GFF3 gene annotations.

The package is designed primarily for genomic analysis in crop species, with soybean examples and Williams 82 genome annotations provided as part of the workflow.

---

## Features

### Global LD decay

Calculates pairwise linkage disequilibrium as squared Pearson correlations (`r²`) between markers within a user-defined physical-distance window.

The fitted LD-decay model is:

$$
r^2 = \frac{1}{1 + Cd}
$$

where $d$ is the physical distance between markers and $C$ is the fitted decay parameter.

### Empirical binned LD decay

Groups pairwise LD observations into user-defined physical-distance intervals and calculates summary statistics such as mean or median LD. This provides a cleaner visualization when millions of pairwise observations are available.

### Regional LD visualization

Generates Haploview-style LD plots and regional LD heatmaps for selected genomic regions or subsets of markers.

### SNP-specific LD analysis

Calculates and visualizes LD surrounding a focal SNP within a user-defined physical window.

### SNP-to-gene mapping

Maps significant or target SNPs to nearby genes using GFF3 genome annotations. The function identifies genes within a user-defined physical window and reports their distance and orientation relative to each SNP.

### LD and gene annotation integration

Combines regional LD patterns with genomic gene models to create publication-quality figures showing LD structure, gene positions, candidate gene annotations, and significant SNP locations.

### SoyBase GFF3 download

Provides a convenient function for downloading and decompressing the Williams 82 soybean GFF3 annotation from SoyBase.

---

# Installation

## Install the development version from GitHub

The development version can be installed directly from GitHub using `devtools`.

```r
# Install devtools if necessary
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install LDdecay
devtools::install_github("shamim-mj/LDdecay")
```
# After installation: 
```r
library(LDdecay)
```
# Input Data
## Genotype data
The genotype dataset should contain:

- One column named taxa identifying samples. Must be the first column!
- Remaining columns containing SNP identifiers.
- Genotype values coded numerically, typically as 0, 1, and 2.

### Genotype Data

The genotype data frame must contain a column named `taxa`, followed by columns containing SNP marker identifiers. Genotype values should be coded as `0`, `1`, or `2`.

| taxa   | ss715000001 | ss715000002 | ss715000003 |
|--------|-------------|-------------|-------------|
| Line001 | 0           | 1           | 2           |
| Line002 | 1           | 1           | 2           |
| Line003 | 0           | 2           | 1           |

### Physical Map

The marker map must contain the following columns:

- `SNP` — variant identifier
- `Chromosome` — chromosome identifier
- `Position` — physical position in base pairs

The structure should resemble:

| SNP         | Chromosome | Position    |
|-------------|------------|-------------|
| ss715000001 | Gm01       | 1562        |
| ss715000002 | Gm01       | 2450        |
| ss715000003 | Gm01       | 1022        |


# GFF3 annotation

GFF3 files should contain standard genomic feature information, including:

- chromosome or sequence name
- start position
- end position
- strand
- feature type
- gene identifiers and annotation attributes

The package primarily uses records where:


```r
type == "gene"
```

# Core Workflows
## 1. Load example data

The package includes example genotype, map, and phenotype datasets.


```r 
library(LDdecay)

data(geno)
data(map)
data(pheno)
```
You can inspect the datasets using:
```r
head(geno)
head(map)
head(pheno)
```

## 2. Calculate Global LD Decay

plot_global_ld_decay() calculates pairwise LD between markers located on the same chromosome and within a specified physical distance.

```r
global_results <- plot_global_ld_decay(
  geno = geno,
  map = map,
  max_distance_bp = 500000,
  max_plot_points = 100000,
  r2_threshold = 0.2,
  output_image_path = "Figures/ld_global_decay_curve.png"
)
```
View the plot:

```r
global_results$plot
```
The complete calculated pairwise LD data are available through:

```r
head(global_results$ld_data)
```
The returned ld_data contains physical distance and corresponding r² values.


## 3. Generate Empirical Binned LD Decay Curves

For large datasets, plotting every pairwise LD observation can produce extremely dense figures.

plot_binned_ld_global_decay() groups observations into physical-distance intervals.

For example, to use 10-kb intervals:

```r
binned_results <- plot_binned_ld_global_decay(
  geno = geno,
  map = map,
  bin_size_bp = 10000,
  max_distance_bp = 500000,
  r2_threshold = 0.2,
  output_image_path = "Figures/ld_binned_decay_profile.png"
)
```
View the resulting plot:

```r
binned_results$plot
```
This approach is useful when millions of pairwise LD observations are generated.

## 4. Create a Haploview-Style LD Plot

plot_haploview_style_ld() creates a regional LD visualization similar to traditional Haploview-style LD matrices.

For example:

```r
plot_haploview_style_ld(
  geno = geno,
  map = map,
  snp_subset = 1:60,
  title = "Soybean Regional LD Matrix"
)
```
The snp_subset argument can be used to specify either marker indices or a selected group of SNP identifiers, depending on the function configuration.


## 5. Calculate LD Around a Focal SNP

plot_single_snp_ld() calculates LD surrounding a focal SNP within a specified physical window.

For example:

```r
single_ld <- plot_single_snp_ld(
  geno = geno,
  map = map,
  focal_snp = "ss715620779",
  window_bp = 100000,
  marker_col = "SNP",
  chr_col = "Chromosome",
  pos_col = "Position",
  min_maf = 0.05
)
```

The returned object can be inspected using:

```r

names(single_ld)
```
For example, the LD matrix can be accessed through:

```r
single_ld$r2
```

## 6. Download the SoyBase Williams 82 GFF3 Annotation

get_soy_gff_from_soybase() downloads the Williams 82 soybean genome annotation from SoyBase, decompresses the file, and imports it into an R data frame.

```r
gff_table <- get_soy_gff_from_soybase(
  dest_dir = "data",
  overwrite = FALSE
)
```
Inspect the annotation:

```r
head(gff_table)
```
Examine available fields:

```r
names(gff_table)
```
The downloaded GFF3 file is retained locally so that it does not need to be downloaded again when overwrite = FALSE.

## 7. Import a Local GFF3 File

If you already have a GFF3 file, use get_gff():

```r

gff_table <- get_gff(
  gff_path = "data/my_annotation.gff3"
)
```
The function converts the imported GFF3 annotation into a standard R data frame.

Inspect the result:

```r
head(gff_table)
```
## 8. Map SNPs to Nearby Genes

sn​p_to_gene_mapping_using_gff3_annot() identifies genes located within a specified physical window around target SNPs.

For example:

```
test_snps <- c(
  "ss715620779",
  "ss715620778",
  "ss715620770"
)
mapped_genes <- snp_to_gene_mapping_using_gff3_annot(
  snp_list = test_snps,
  gff_path = "data/glyma.Wm82.gnm1.ann1.DvBy.gene_models_main.gff3",
  map_data = map,
  window_bp = 50000
)
```
Inspect the results:

```r
head(mapped_genes)
```

The returned table includes information such as:

- query SNP
- chromosome
- SNP position
- gene ID
- gene symbol
- gene coordinates
- strand
- distance from the SNP
- upstream/downstream/overlapping orientation
- gene description

For example:

query_snp
snp_chr
snp_pos
gene_id
gene_symbol
start
end
strand
distance_to_snp
orientation
description


## 9. Plot LD and Gene Annotations in a Regional Area of Interest

plot_snp_ld_region_and_genes() combines regional LD information with gene models from a GFF3 annotation.

This produces a composite figure containing:

Gene models across the genomic region.
Candidate genes matching a user-defined annotation pattern.
Significant SNP locations.
Regional pairwise LD.
Physical genomic coordinates.

Example:

```r
result <- plot_snp_ld_region_and_genes(
  map = map,
  geno = geno,
  gff_table = gff_table,
  chromosome = "Gm03",
  focal_snp = "ss715620779",
  ld_window = 50000,
  gene_name = "ABC/PDR Transport",
  gene_pattern = "ABC-2/plant PDR ABC transporter|ABC|PDR transporter"
)
```

The resulting object contains several useful components:

```r
result$plot
result$gene_plot
result$ld_plot
result$genes
result$cyp_genes
result$marker_map
result$significant_map
result$r2
result$ld_data
```
Display the complete figure: It will be automatically displayed without even running this code

```r
result$plot
```
# Candidate Gene Filtering

The regional plotting function allows users to search gene annotations using regular expressions.

For example, to identify cytochrome P450-related annotations:

```r
gene_pattern <- "cytochrome p450|cytochrome|p450"
```
Then:

```r
result <- plot_snp_ld_region_and_genes(
  map = map,
  geno = geno,
  gff_table = gff_table,
  chromosome = "Gm03",
  focal_snp = "ss715620779",
  ld_window = 50000,
  gene_pattern = gene_pattern,
  gene_name = "Cytochrome P450"
)
```
This allows the same visualization framework to be used for different biological pathways or candidate gene families.

# Significant SNPs

A vector of significant SNPs can be supplied using significant_snps.

For example:

```r
significant_snps <- c(
  "ss715620779",
  "ss715620778",
  "ss715620770"
)
```
then

```r
result <- plot_snp_ld_region_and_genes(
  map = map,
  geno = geno,
  gff_table = gff_table,
  chromosome = "Gm03",
  focal_snp = "ss715620779",
  ld_window = 50000,
  significant_snps = significant_snps
)
```
Significant SNPs are shown as vertical reference lines in the regional plot.

# Choosing the LD Window

The physical window can be adjusted depending on the objective.

For example:

```r

ld_window = 50000
```
uses a 50-kb window around the focal SNP.

A larger window can be used for broader regional exploration:

```r
ld_window = 360000
```
This corresponds to a 360-kb flanking window.

The appropriate window depends on the species, population, marker density, and expected LD decay.


# Output Files

Most plotting functions allow the user to save figures directly.

For example:

```r
plot_global_ld_decay(
  geno = geno,
  map = map,
  output_image_path = "Figures/global_ld_decay.png"
)
```
The output directory is automatically created if it does not already exist.

Figures are saved at publication-quality resolution.

# Package Documentation

Full function documentation is available through the package website:

[LDdecay Documentation](https://shamim-mj.github.io/LDdecay/)

The documentation includes:

function descriptions
argument definitions
return values
examples
package workflows

# Development

The source code and development versions of LDdecay are available on GitHub:

https://github.com/shamim-mj/LDdecay

Bug reports and suggestions are welcome through the GitHub repository.

# License

LDdecay is distributed under the MIT License. See the [LICENSE](LICENSE) file for details.

# Citation

If you use LDdecay in research or publications, please cite the package according to the citation information provided in the repository and package documentation.

Run:

```r
citation(LDdecay)
```







