\# LDdecay: Advanced Tools for Linkage Disequilibrium Decay and Functional Mapping



`LDdecay` is an optimized R package designed for processing, analyzing, and generating publication-ready visualizations for global genome-wide Linkage Disequilibrium (LD) decay, binned empirical decay curves, regional LD block structure maps, and fast overlapping variant-to-gene mapping from standardized GFF3 annotations.



\## Features

\- \*\*Global LD Decay Model\*\*: Fits the classical genetic non-linear least squares (NLS) regression curve ($r^2 \\sim 1/(1 + C \\times d)$) over millions of marker combinations using an optimized sliding physical window.

\- \*\*Empirical Binned Decay Curves\*\*: Smooths and groups dense pairwise data into custom intervals (e.g., 10kb windows) using Cartesian coordinate zooms to entirely avoid edge boundaries or cropping spikes.

\- \*\*Regional Block Plots \& Heatmaps\*\*: Generates publication-ready Haploview-style inverted triangle matrices and local area-of-interest annotations flanking candidate regions.

\- \*\*Ultra-Fast Variant Feature Scanner\*\*: Maps target variant coordinates to GFF3 structural genomic annotations using vector-optimized memory management.



\---



\## Installation



You can install the development version of `LDdecay` straight from GitHub using `devtools`:



```R

\# Install devtools if you haven't already

if (!requireNamespace("devtools", quietly = TRUE)) {

&#x20; install.packages("devtools")

}



\# Install LDdecay directly from your public GitHub repository

devtools::install\_github("yourusername/LDdecay")

```



\---



\## Core Package Workflows \& Examples



Below are standalone script examples showing how to run the package functions. The package includes built-in test data arrays (`geno`, `map`, `pheno`, and a subsetted chromosome 3 GFF3 file) to let you run tests immediately.



\### 1. Load Package Data

```R

library(LDdecay)

library(ggplot2)



\# Load package internal data subsets 

data(geno)  # Row coordinates are taxa samples, columns are 0/1/2 SNP values

data(map)   # Variant mapping coordinates (SNP, Chromosome, Position)

data(pheno) # Phenotype tracking matrices

```



\### 2. Calculate and Plot Global LD Decay Curve (NLS Regression Line)

This function isolates short-range regions, downsamples massive matrices safely to prevent system RAM crashes, fits an NLS model, and renders the scatter decay trend line.



```R

\# Calculate global linkage profiles and save plot asset to file

global\_results <- plot\_global\_ld\_decay(

&#x20; geno = geno, 

&#x20; map = map, 

&#x20; max\_distance\_bp = 500000,    # Focus window on 500kb tail

&#x20; max\_plot\_points = 100000,    # Max background scatter points

&#x20; r2\_threshold = 0.2,          # Baseline threshold marker line

&#x20; output\_image\_path = "Figures/ld\_global\_decay\_curve.png"

)



\# View the ggplot object directly in RStudio

print(global\_results\\$plot)

```



\### 3. Generate Empirical Binned Decay Curves

Groups raw data combinations into fixed base-pair bins to calculate mean/median trends without crowding plots with millions of overlapping data points.



```R

binned\_results <- plot\_binned\_ld\_decay(

&#x20; geno = geno, 

&#x20; map = map, 

&#x20; bin\_size\_bp = 10000,         # 10kb physical grouping intervals

&#x20; max\_distance\_bp = 500000,    # Focus up to 500kb

&#x20; r2\_threshold = 0.2,

&#x20; output\_image\_path = "Figures/ld\_binned\_decay\_profile.png"

)



\# View summary trend lines

print(binned\_results\\$plot)

```



\### 4. Create a Haploview-Style Triangle Heatmap

Subsets your target physical coordinate structures into a clean square correlation matrix array and prints traditional genetic heat blocks.



```R

\# Extract an automated slice of the first 60 rows, or declare index vectors

plot\_haploview\_ld(

&#x20; geno = geno, 

&#x20; map = map, 

&#x20; snp\_subset = 1:60,           # Accept row indices (e.g., 1:60) or specific marker vectors

&#x20; title = "Soybean Regional Block Matrix"

)

```



\### 5. High-Speed Variant Feature Extraction mapping via GFF3

Pulls structural coordinates from `.gff3` layers and returns clean upstream, downstream, or overlapping orientations alongside physical genomic annotations.



```R

\# Locate the internal GFF3 text file pathway dynamically

gff\_file\_path <- system.file("extdata", "glyma\_soy\_models\_chr3.gff3", package = "LDdecay")



\# Extract overlapping gene attributes within 50kb physical windows flanking target SNPs

mapped\_genes\_table <- get\_bulk\_snp\_genes(

&#x20; snp\_list = c("ss715585023", "ss715629879"), 

&#x20; gff\_path = gff\_file\_path, 

&#x20; map\_data = map, 

&#x20; window\_bp = 50000

)



\# Preview polished dataset results

head(mapped\_genes\_table)

```



\---



\## Technical Input Guidelines



To ensure your datasets run seamlessly through the package pipelines, format your files as follows:



\- \*\*Genotype (`geno`)\*\*: Must include a column named exactly `"taxa"` containing unique string sample labels. All other columns must match your structural SNP variant marker keys exactly, scaled as numeric values or character codes (`0`, `1`, `2`).

\- \*\*Physical Map (`map`)\*\*: A standard data frame requiring three explicit column titles case-sensitive: `"SNP"`, `"Chromosome"`, and `"Position"` (numeric base pairs).

\- \*\*GFF3 (`.gff3`)\*\*: Standard feature structural file schemas containing records where `type == "gene"` populated with attribute bindings for `ID`, `Name`, and `Note`.



\## License

Distributed under the MIT License. See `LICENSE` for more information.



