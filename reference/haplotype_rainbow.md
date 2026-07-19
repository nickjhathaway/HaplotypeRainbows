# Create a HaplotypeRainbow

Convenience wrapper around `HaplotypeRainbow$new()`.

## Usage

``` r
haplotype_rainbow(
  data,
  sample_col = "library_sample_name",
  target_col = "target_name",
  popuid_col = "seq",
  rel_abund_col = "reads"
)
```

## Arguments

- data:

  A data frame / tibble with one row per (sample, target, haplotype).

- sample_col:

  Name of the sample identifier column.

- target_col:

  Name of the target / locus column.

- popuid_col:

  Name of the within-target haplotype identifier column.

- rel_abund_col:

  Name of the relative counts column.

## Value

A
[HaplotypeRainbow](https://nickjhathaway.github.io/HaplotypeRainbows/reference/HaplotypeRainbow.md)
object.

## Examples

``` r
if (FALSE) { # \dontrun{
rb <- haplotype_rainbow(my_allele_table)
} # }
```
