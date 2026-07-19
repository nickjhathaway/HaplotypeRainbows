# Example target metadata for [pfIsosHeomeV1](https://nickjhathaway.github.io/HaplotypeRainbows/reference/pfIsosHeomeV1.md)

Per-target metadata for the targets in
[pfIsosHeomeV1](https://nickjhathaway.github.io/HaplotypeRainbows/reference/pfIsosHeomeV1.md),
for demonstrating the target-annotation strip. The `target_name` column
matches `pfIsosHeomeV1`'s `p_name`.

## Usage

``` r
pfIsosHeomeV1_targetMeta
```

## Format

A data frame with 100 rows (one per target):

- target_name:

  Target / locus name (matches `p_name`).

- gene_id:

  Overlapping gene identifier.

- gene_description:

  Gene description.

- group:

  Whether the target is included or excluded ("include"/"exclude").

- class:

  Target class ("Diversity"/"Drug").

- chrom:

  Chromosome.
