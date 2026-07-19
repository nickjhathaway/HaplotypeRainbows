# Example isolate haplotype table with alternative column names

The same data as
[pfIsosHeomeV1](https://nickjhathaway.github.io/HaplotypeRainbows/reference/pfIsosHeomeV1.md)
but with alternative column names, used to demonstrate supplying a
custom column mapping to
[HaplotypeRainbow](https://nickjhathaway.github.io/HaplotypeRainbows/reference/HaplotypeRainbow.md).

## Usage

``` r
pfIsosHeomeV1_otherName
```

## Format

A tibble with 7,611 rows and 4 columns:

- Sample:

  Sample identifier.

- loci:

  Target / locus name.

- ID:

  Within-target haplotype identifier.

- freq:

  Within-sample relative abundance of the haplotype.

## Source

SeekDeep output; see <https://github.com/bailey-lab/SeekDeep>.
