#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr %>%
#' @importFrom rlang := sym
#' @importFrom R6 R6Class
#' @importFrom stats hclust dist cutree as.dendrogram
#' @importFrom tidyr pivot_wider
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggnewscale new_scale_fill
#' @importFrom grDevices colorRampPalette cairo_pdf pdf dev.off col2rgb rgb2hsv
#' @importFrom grid grid.newpage grid.draw
#' @importFrom gtable gtable_filter
## usethis namespace: end
NULL

# The prep engine reuses the original dplyr pipelines, which reference many derived
# columns by bare name via non-standard evaluation. Declare them so R CMD check does
# not flag them as undefined global variables.
utils::globalVariables(c(
  ".data",
  # canonical column names used inside the engine
  "sample", "target", "hapid", "rel_abund",
  # derived normalisation / geometry columns
  "n_targets", "total_abund", "within_sample_freq", "within_sample_freq_mod",
  "sample_coi", "freq_cumsum", "freq_mod_cumsum",
  "fake_freq", "fake_freq_mod", "fake_freq_cumsum", "fake_freq_mod_cumsum",
  # ranking columns
  "samp_n", "pop_id", "max_pop_id", "pop_id_frac", "hue_mod",
  "pop_id_perc", "pop_id_frac_reg_color", "pop_id_perc_log",
  "pop_id_frac_log_color", "is_invariant",
  # shade colouring columns
  "n", "total", "freq", "p_uniq_haps", "h_id", "h_id_freq", "h_id_freq_mod",
  "cum_freq", "mod_cum_freq", "p_color_id", "p_hue", "p_base_color",
  "h_color", "h_color_mod", "h_color_by_freq", "h_color_by_freq_mod",
  # clustering / reorder helpers
  "new_identifier", "marker", "sample_count", "sample_freq"
))
