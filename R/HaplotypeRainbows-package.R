#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr %>%
#' @importFrom rlang := sym
#' @importFrom R6 R6Class
#' @importFrom stats hclust dist
#' @importFrom tidyr pivot_wider
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggnewscale new_scale_fill
#' @importFrom grDevices colorRampPalette cairo_pdf pdf dev.off col2rgb rgb2hsv
## usethis namespace: end
NULL

# The prep engine reuses the original dplyr pipelines, which reference many derived
# columns by bare name via non-standard evaluation. Declare them so R CMD check does
# not flag them as undefined global variables.
utils::globalVariables(c(
  ".data",
  # canonical column names used inside the engine
  "s_Sample", "p_name", "h_popUID", "c_AveragedFrac",
  # derived geometry / ranking columns
  "targetNumber", "totalAbund", "s_COI", "relAbundCol_mod",
  "fracCumSum", "fracModCumSum", "fakeFrac", "fakeFracMod",
  "fakeFracCumSum", "fakeFracModCumSum", "samp_n", "popid", "maxPopid",
  "popidFrac", "hueMod", "popidPerc", "popidFracRegColor",
  "popidPercLog", "popidFracLogColor", "is_invariant",
  # shade colouring columns
  "n", "total", "freq", "p_uniqHaps", "h_id", "h_id_freq", "h_id_freq_mod",
  "cumFreq", "modCumFreq", "p_color_ID", "p_hue", "p_baseColor",
  "h_color", "h_color_mod", "h_color_byFreq", "h_color_byFreq_mod",
  # clustering / reorder helpers
  "new_identifier", "marker", "sample_count", "sample_freq"
))
