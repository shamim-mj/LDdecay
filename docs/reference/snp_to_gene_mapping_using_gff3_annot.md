# Map SNPs to Nearby Genes Using GFF3 Genome Annotations

Maps a set of target SNPs to nearby gene models using genomic
coordinates from a GFF3 annotation file. For each SNP, the function
searches upstream and downstream within a user-defined physical window
and identifies all genes whose genomic intervals overlap that window.

## Usage

``` r
snp_to_gene_mapping_using_gff3_annot(
  snp_list,
  gff_path,
  map_data,
  window_bp = 50000
)
```

## Arguments

- snp_list:

  Character vector of target SNP identifiers to query. Typically, this
  is a vector containing significant SNPs identified from a GWAS or
  other marker-trait association analysis.

- gff_path:

  Character string specifying the path to a local GFF3 annotation file.
  The file should contain genomic coordinates and gene annotations. For
  the soybean Williams 82 reference genome, annotations can be obtained
  from the SoyBase genome annotation repository:
  <https://data.soybase.org/Glycine/max/annotations/Wm82.gnm1.ann1.DvBy/>.

- map_data:

  A data frame containing physical marker coordinates. It must contain
  columns named `SNP`, `Chromosome`, and `Position`. `SNP` must contain
  marker identifiers matching those supplied in `snp_list`.

- window_bp:

  Numeric value specifying the upstream and downstream physical search
  distance around each SNP, in base pairs. A gene is retained when its
  genomic interval overlaps the SNP-centered window. Default is `50000`
  (50 kb).

## Value

A tibble containing one row for each SNP-gene association found within
the specified physical window. The returned table contains:

- `query_snp`: Target SNP identifier.

- `snp_chr`: Chromosome containing the SNP.

- `snp_pos`: Physical position of the SNP in base pairs.

- `gene_id`: Gene identifier extracted from the GFF3 annotation.

- `gene_symbol`: Gene name or symbol extracted from the GFF3 `Name`
  attribute, when available.

- `seqnames`: Original chromosome/sequence identifier from the GFF3
  file.

- `start`: Gene start coordinate.

- `end`: Gene end coordinate.

- `strand`: Gene strand orientation.

- `type`: GFF3 feature type; retained as `"gene"`.

- `distance_to_snp`: Distance in base pairs between the SNP position and
  the midpoint of the gene. Positive values indicate that the gene
  midpoint is downstream of the SNP, while negative values indicate that
  it is upstream.

- `orientation`: Relative position of the gene with respect to the SNP:
  `"Upstream"`, `"Downstream"`, or `"Overlapping"`.

- `description`: Functional or descriptive annotation extracted from the
  GFF3 `Note` field, when available.

## Details

The function is designed for post-GWAS candidate-gene analysis. It first
synchronizes SNP identifiers and genomic coordinates with the supplied
physical map, imports the GFF3 annotation using rtracklayer, and
extracts only features classified as genes. Gene coordinates and
selected annotation metadata are then used to identify genes located
within the specified distance of each target SNP.

Gene position relative to each SNP is classified as `"Upstream"`,
`"Downstream"`, or `"Overlapping"` based on the position of the gene
midpoint relative to the SNP coordinate.

This function currently expects chromosome identifiers that contain the
soybean chromosome naming convention `"Gm01"`–`"Gm20"`. It is therefore
primarily intended for soybean GFF3 annotations unless the
chromosome-extraction step is modified for another genome.

The mapping procedure consists of four main steps:

1.  Target SNP identifiers are cleaned and matched against the supplied
    physical map.

2.  The GFF3 file is imported and filtered to retain only features with
    `type == "gene"`. Selected attributes, including gene ID, gene name,
    and description, are extracted into a lightweight data frame.

3.  SNPs and genes are matched by chromosome, followed by an interval
    overlap test. A gene is retained when its start coordinate is before
    the end of the SNP-centered search window and its end coordinate is
    after the beginning of that window.

4.  The distance between each SNP and the midpoint of each nearby gene
    is calculated and used to classify the gene as upstream, downstream,
    or overlapping.

The function does not modify the original genotype or map data. It
returns only SNP-gene relationships detected within the specified
physical window.

## Interpretation

A gene is considered "Overlapping" when its midpoint is exactly equal to
the SNP position. This classification is based on gene midpoint rather
than gene interval overlap. A gene can therefore physically overlap the
SNP while still being classified as upstream or downstream if its
midpoint lies on one side of the SNP.

## Genome compatibility

Chromosome identifiers are standardized by extracting patterns matching
`Gm\d+`. Consequently, this implementation is optimized for soybean
chromosome naming conventions such as `Gm01`, `Gm02`, and `Gm20`. For
other crops, the chromosome-standardization step should be adapted to
the naming convention used by the corresponding GFF3 file.

## Examples

``` r
if (FALSE) { # \dontrun{
library(LDdecay)

# Significant SNPs identified from a GWAS
test_snps <- c(
  "ss715620779",
  "ss715620778",
  "ss715620770"
)

# Map SNPs to genes within 50 kb
test_mapping <- snp_to_gene_mapping_using_gff3_annot(
  snp_list = test_snps,
  gff_path = "data/glyma.Wm82.gnm1.ann1.DvBy.gene_models_main.gff3",
  map_data = map,
  window_bp = 50000
)

# Inspect the SNP-gene associations
head(test_mapping)

# View genes surrounding each target SNP
print(test_mapping)
} # }
```
