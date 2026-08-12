# Download the SoyBase Williams 82 GFF3 Annotation

Downloads the Williams 82 soybean genome annotation from SoyBase,
decompresses the GFF3 archive, and converts the annotation to a standard
data frame using [`get_gff`](get_gff.md).

## Usage

``` r
get_soy_gff_from_soybase(dest_dir = "data", overwrite = FALSE)
```

## Arguments

- dest_dir:

  Character string specifying the directory in which the GFF3 annotation
  file will be downloaded and stored. The directory is created
  recursively if it does not already exist. Default is `"data"`.

- overwrite:

  Logical indicating whether an existing uncompressed GFF3 file should
  be replaced. Default is `FALSE`.

## Value

A `data.frame` containing the parsed Williams 82 soybean genome
annotation.

## Details

This function is specifically designed for soybean (*Glycine max*) and
retrieves the Williams 82 genome annotation from SoyBase.

If the uncompressed GFF3 file already exists in `dest_dir` and
`overwrite = FALSE`, the existing file is used and no download is
performed.

The function downloads
`glyma.Wm82.gnm1.ann1.DvBy.gene_models_main.gff3.gz` from SoyBase,
decompresses the archive, removes the compressed file, and passes the
resulting GFF3 file to [`get_gff`](get_gff.md) for parsing.

## References

SoyBase. Soybean Genome Database. <https://www.soybase.org/>

## Examples

``` r
if (FALSE) { # \dontrun{
library(LDdecay)

gff_table <- get_soy_gff_from_soybase(
  dest_dir = "data"
)

head(gff_table)
} # }
```
