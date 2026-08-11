# LDdecay: Advanced Tools for Linkage Disequilibrium Decay and Functional Mapping
# LDdecay <img src="man/figures/logo.png" align="right" height="139" />

<!-- Badges start here -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://repostatus.org)](https://repostatus.org)
[![Lifecycle: Experimental](https://shields.io)](https://r-lib.org)
[![Minimal R Version](https://shields.io)](https://r-project.org)
<!-- Badges end here -->



`LDdecay` is an optimized R package designed for processing, analyzing, and generating publication-ready visualizations for global genome-wide Linkage Disequilibrium (LD) decay, binned empirical decay curves, regional LD block structure maps, and fast overlapping variant-to-gene mapping from standardized GFF3 annotations.



## Features

 **Global LD Decay Model**: Fits the classical genetic non-linear least squares (NLS) regression curve ($r^2 \\sim 1/(1 + C \\times d)$) over millions of marker combinations using an optimized sliding physical window.

**Empirical Binned Decay Curves**: Smooths and groups dense pairwise data into custom intervals (e.g., 10kb windows) using Cartesian coordinate zooms to entirely avoid edge boundaries or cropping spikes.

**Regional Block Plots & Heatmaps**: Generates publication-ready Haploview-style inverted triangle matrices and local area-of-interest annotations flanking candidate regions.

**Ultra-Fast Variant Feature Scanner**: Maps target variant coordinates to GFF3 structural genomic annotations using vector-optimized memory management.







## Installation



You can install the development version of `LDdecay` straight from GitHub using `devtools`:



```R

# Install devtools if you haven't already

if (!requireNamespace("devtools", quietly = TRUE)) {

install.packages("devtools")

}



# Install LDdecay directly from your public GitHub repository

devtools::install_github("shamim-mj/LDdecay")

```







## Core Package Workflows & Examples



Below are standalone script examples showing how to run the package functions. The package includes built-in test data arrays (`geno`, `map`, `pheno`, and a subsetted chromosome 3 GFF3 file) to let you run tests immediately.



### 1. Load Package Data

```R

library(LDdecay)

library(ggplot2)



# Load package internal data subsets 

data(geno)  # Row coordinates are taxa samples, columns are 0/1/2 SNP values

data(map)   # Variant mapping coordinates (SNP, Chromosome, Position)

data(pheno) # Phenotype tracking matrices

```



### 2. Calculate and Plot Global LD Decay Curve (NLS Regression Line)

This function isolates short-range regions, downsamples massive matrices safely to prevent system RAM crashes, fits an NLS model, and renders the scatter decay trend line.



```R

# Calculate global linkage profiles and save plot asset to file

global_results <- plot_global_ld_decay(

geno = geno, 

map = map, 

max_distance_bp = 500000,    # Focus window on 500kb tail

max_plot_points = 100000,    # Max background scatter points

r2_threshold = 0.2,          # Baseline threshold marker line

output_image_path = "Figures/ld_global_decay_curve.png"

)



# View the ggplot object directly in RStudio

print(global_results$plot)

```



### 3. Generate Empirical Binned Decay Curves

Groups raw data combinations into fixed base-pair bins to calculate mean/median trends without crowding plots with millions of overlapping data points.



```R

binned_results <- plot_binned_ld_global_decay(

geno = geno, 

map = map, 

bin_size_bp = 10000,         # 10kb physical grouping intervals

max_distance_bp = 500000,    # Focus up to 500kb

r2_threshold = 0.2,

output_image_path = "Figures/ld_binned_decay_profile.png"

)



# View summary trend lines

print(binned_results$plot)

```



### 4. Create a Haploview-Style Triangle Heatmap

Subsets your target physical coordinate structures into a clean square correlation matrix array and prints traditional genetic heat blocks.



```R

# Extract an automated slice of the first 60 rows, or declare index vectors

plot_haploview_style_ld(

geno = geno, 

map = map, 

snp_subset = 1:60,           # Accept row indices (e.g., 1:60) or specific marker vectors

title = "Soybean Regional Block Matrix"

)

```


### 5. Plot Linkage Disequilibrium Heatmap Around a Focal Variant

Isolates flanking markers inside a physical window, filters out variants with a low Minor Allele Frequency (MAF), and renders a correlation matrix heatmap.



```R

# A list containing ggplot object, the raw R2 correlation matrix, and filtered maps/genotypes.
# This plot shows the correlation of the significant or focal SNP to other SNPs. Focal SNP is indicated with a "*"

plot_single_snp_ld(

geno = geno, 
map = map,
window_bp = 100000,
marker_col = "SNP",
chr_col = "Chromosome",
pos_col = "Position"
 min_maf = 0.05
snp_subset = 1:60)

```


### 6. Download, Decompress, and Import SoyBase GFF3 Annotation Files

Automatically connects to the authoritative SoyBase repository data collection, downloads the compressed Williams 82 genome annotation layer (.gff3.gz), and 
extracts it into a memory-optimized data frame ready for variant structural mapping.


```R

# Give the destination folder a name!


# Extract overlapping gene attributes within 50kb physical windows flanking target SNPs

gff_table <- get_soy_gff_from_soybase(

dest_dir = "data",
overwrite = FALSE
)



# Preview polished dataset results

head(gff_table)

```






### 7. High-Speed Variant Feature Extraction mapping via GFF3

Pulls structural coordinates from `.gff3` layers and returns clean upstream, downstream, or overlapping orientations alongside physical genomic annotations.



```R

# Locate the internal GFF3 text file pathway dynamically


# Extract overlapping gene attributes within 50kb physical windows flanking target SNPs

mapped_genes_table <- snp_to_gene_mapping_using_gff3_annot(

snp_list = c("ss715585023", "ss715629879"), 

gff_path = gff_file_path, 

map_data = map, 

window_bp = 50000

)



# Preview polished dataset results

head(mapped_genes_table)

```


### 8. Plot Linkage Disequilibrium Matrix and Gene Annotations in an Area of Interest (AOI)

Renders a publication-quality composite plot aligning a pairwise Linkage Disequilibrium (LD) correlation heatmap with structural gene models extracted from a GFF database.
Highlights candidate gene sets matching key functional patterns and tracks user-defined variant subsets.



```R

# Locate the GFF3 table, map file, and create a victor of your significant SNPs


# Extract overlapping gene attributes within 50kb physical windows flanking target SNPs

plot_snp_ld_region_and_genes(

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
)


```








## Technical Input Guidelines



To ensure your datasets run seamlessly through the package pipelines, format your files as follows:



**Genotype (`geno`)**: Must include a column named exactly `"taxa"` containing unique string sample labels. All other columns must match your structural SNP variant marker keys exactly, scaled as numeric values or character codes (`0`, `1`, `2`).

**Physical Map (`map`)**: A standard data frame requiring three explicit column titles case-sensitive: `"SNP"`, `"Chromosome"`, and `"Position"` (numeric base pairs).

**GFF3 (`.gff3`)**: Standard feature structural file schemas containing records where `type == "gene"` populated with attribute bindings for `ID`, `Name`, and `Note`.



## License

Distributed under the MIT License. See `LICENSE` for more information.



