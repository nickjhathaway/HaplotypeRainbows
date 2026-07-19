# HaplotypeRainbows: Creates Haplotype Rainbow Plots

Creates haplotype "rainbow" plots from targeted amplicon
(microhaplotype) sequencing data. Given a table of samples, targets,
within-target haplotype identifiers and relative counts, it builds a
grid where columns are targets and rows are samples, with each cell
showing the within-sample haplotype composition. Haplotype colours
rotate across targets to produce the characteristic rainbow. The package
exposes a single R6 class, HaplotypeRainbow, that carries the column
mapping so it only has to be set once. Column defaults follow the
Portable Microhaplotype Object (PMO) convention.

## See also

Useful links:

- <https://github.com/nickjhathaway/HaplotypeRainbows>

- <https://nickjhathaway.github.io/HaplotypeRainbows/>

- Report bugs at
  <https://github.com/nickjhathaway/HaplotypeRainbows/issues>

## Author

**Maintainer**: Nicholas Hathaway <nickjhathaway@gmail.com>
([ORCID](https://orcid.org/0000-0001-9639-2894))

Authors:

- Nicholas Hathaway <nickjhathaway@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-9639-2894))
