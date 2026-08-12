# Import and Convert GFF3 Annotation Files to Standard Data Frames

Reads a local genomic annotation file in GFF3 format and converts it
into a standard R data frame for downstream genomic analyses.

## Usage

``` r
get_gff(gff_path)
```

## Arguments

- gff_path:

  Character string specifying the path to a GFF3 annotation file.

## Value

A `data.frame` containing parsed genomic features, including chromosome,
start and end coordinates, strand, feature type, and annotation
attributes.

## Examples

``` r
if (FALSE) { # \dontrun{
library(LDdecay)

gff_table <- get_gff("data/annotation.gff3")
head(gff_table)
} # }
```
