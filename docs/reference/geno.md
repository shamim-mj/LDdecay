# Example Genotype Dataset

A small soybean genotype dataset included with LDdecay for demonstrating
linkage disequilibrium analysis and visualization.

## Usage

``` r
geno
```

## Format

A data frame containing example soybean lines and SNP markers. The
`taxa` column identifies the genotype, and the remaining columns contain
SNP genotype dosages coded as 0, 1, or 2.

## Examples

``` r
data(geno)
head(geno)
#>                taxa ss715620471 ss715620532 ss715620606 ss715620667 ss715620670
#> PI417398   PI417398           0           0           0           0           0
#> PI91750     PI91750           0           0           0           0           0
#> PI548614   PI548614           0           0           0           0           0
#> PI594457B PI594457B           2           2           2           0           0
#> PI406684   PI406684           2           2           2           0           0
#> PI423973   PI423973           0           0           0           0           0
#>           ss715620714 ss715620779 ss715620815 ss715620829 ss715620870
#> PI417398            0           0           0           0           0
#> PI91750             0           0           2           2           2
#> PI548614            0           2           0           0           0
#> PI594457B           2           2           2           0           0
#> PI406684            0           0           0           0           0
#> PI423973            0           0           0           0           0
#>           ss715620924 ss715620928 ss715620957 ss715621026 ss715621053
#> PI417398            2           2           2           0           0
#> PI91750             0           0           0           0           0
#> PI548614            0           2           0           0           0
#> PI594457B           2           2           2           0           0
#> PI406684            2           2           0           0           0
#> PI423973            2           2           2           0           0
#>           ss715621062 ss715621069 ss715621070
#> PI417398            2           0           0
#> PI91750             0           0           0
#> PI548614            0           0           2
#> PI594457B           2           0           0
#> PI406684            2           0           2
#> PI423973            2           0           0
```
