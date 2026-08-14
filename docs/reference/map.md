# Example Physical Marker Map

A small soybean physical marker map corresponding to the SNP markers
included in the `geno` example dataset.

## Usage

``` r
map
```

## Format

A data frame with three columns:

- SNP:

  SNP marker identifier.

- Chromosome:

  Chromosome identifier.

- Position:

  Physical position in base pairs.

## Examples

``` r
data(map)
head(map)
#>           SNP Chromosome Position
#> 1 ss715620471       Gm15  1270358
#> 2 ss715620532       Gm15  1304136
#> 3 ss715620606       Gm15  1349135
#> 4 ss715620667       Gm15  1387913
#> 5 ss715620670       Gm15  1389521
#> 6 ss715620714       Gm15  1416378
```
