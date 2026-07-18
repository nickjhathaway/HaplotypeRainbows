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
