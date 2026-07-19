# Example *Plasmodium falciparum* isolate haplotype table (SeekDeep column names)

Microhaplotype calls from *Plasmodium falciparum* lab isolates and
mixtures, one row per (sample, target, haplotype). Uses the historical
SeekDeep column names; pass these names explicitly to
[HaplotypeRainbow](https://nickjhathaway.github.io/HaplotypeRainbows/reference/HaplotypeRainbow.md)
since the class now defaults to the PMO convention.

## Usage

``` r
pfIsosHeomeV1
```

## Format

A tibble with 7,611 rows and 4 columns:

- s_Sample:

  Sample identifier.

- p_name:

  Target / locus name.

- h_popUID:

  Within-target haplotype (population) identifier.

- c_AveragedFrac:

  Within-sample relative abundance of the haplotype.

## Source

SeekDeep output; see <https://github.com/bailey-lab/SeekDeep>.
