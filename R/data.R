#' Example *Plasmodium falciparum* isolate haplotype table (SeekDeep column names)
#'
#' Microhaplotype calls from *Plasmodium falciparum* lab isolates and mixtures,
#' one row per (sample, target, haplotype). Uses the historical SeekDeep column
#' names; pass these names explicitly to [HaplotypeRainbow] since the class now
#' defaults to the PMO convention.
#'
#' @format A tibble with 7,611 rows and 4 columns:
#' \describe{
#'   \item{s_Sample}{Sample identifier.}
#'   \item{p_name}{Target / locus name.}
#'   \item{h_popUID}{Within-target haplotype (population) identifier.}
#'   \item{c_AveragedFrac}{Within-sample relative abundance of the haplotype.}
#' }
#' @source SeekDeep output; see <https://github.com/bailey-lab/SeekDeep>.
"pfIsosHeomeV1"

#' Example isolate haplotype table with alternative column names
#'
#' The same data as [pfIsosHeomeV1] but with alternative column names, used to
#' demonstrate supplying a custom column mapping to [HaplotypeRainbow].
#'
#' @format A tibble with 7,611 rows and 4 columns:
#' \describe{
#'   \item{Sample}{Sample identifier.}
#'   \item{loci}{Target / locus name.}
#'   \item{ID}{Within-target haplotype identifier.}
#'   \item{freq}{Within-sample relative abundance of the haplotype.}
#' }
#' @source SeekDeep output; see <https://github.com/bailey-lab/SeekDeep>.
"pfIsosHeomeV1_otherName"

#' Example sample metadata for [pfIsosHeomeV1]
#'
#' Per-sample metadata for the samples in [pfIsosHeomeV1], for demonstrating the
#' sample-metadata sidebar. Subset to the example samples with columns that are
#' constant across the subset dropped. The `sample` column matches `pfIsosHeomeV1`'s
#' `s_Sample`.
#'
#' @format A data frame with 61 rows (one per sample):
#' \describe{
#'   \item{sample}{Sample identifier (matches `s_Sample`).}
#'   \item{collection_date}{Collection date, where known.}
#'   \item{country}{Country of origin.}
#'   \item{site}{Collection site.}
#'   \item{collection_year}{Collection year, where known.}
#'   \item{region}{Geographic region.}
#'   \item{subRegion}{Geographic sub-region.}
#'   \item{secondaryRegion}{Broader (continent-level) region.}
#' }
"pfIsosHeomeV1_sampleMeta"

#' Example target metadata for [pfIsosHeomeV1]
#'
#' Per-target metadata for the targets in [pfIsosHeomeV1], for demonstrating the
#' target-annotation strip. The `target_name` column matches `pfIsosHeomeV1`'s `p_name`.
#'
#' @format A data frame with 100 rows (one per target):
#' \describe{
#'   \item{target_name}{Target / locus name (matches `p_name`).}
#'   \item{gene_id}{Overlapping gene identifier.}
#'   \item{gene_description}{Gene description.}
#'   \item{group}{Whether the target is included or excluded ("include"/"exclude").}
#'   \item{class}{Target class ("Diversity"/"Drug").}
#'   \item{chrom}{Chromosome.}
#' }
"pfIsosHeomeV1_targetMeta"
