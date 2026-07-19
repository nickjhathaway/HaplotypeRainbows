#' Colour-blind-friendly categorical palettes
#'
#' Three categorical palettes (8, 12 and 15 colours) drawn from colour-blind-friendly
#' colour sets, suitable for colouring sample/target metadata bands or for use as the
#' `colors` ramp in [HaplotypeRainbow]'s `plot()`.
#'
#' @format Character vectors of hex colour strings:
#' \describe{
#'   \item{colorPalette_08}{8 colours.}
#'   \item{colorPalette_12}{12 colours.}
#'   \item{colorPalette_15}{15 colours.}
#' }
#' @name color_palettes
#' @examples
#' scales::show_col(colorPalette_12)
NULL

#' @rdname color_palettes
#' @export
colorPalette_08 <- c(
  "#2271B2", "#3DB7E9", "#F748A5", "#359B73",
  "#D55E00", "#E69F00", "#F0E442", "#000000"
)

#' @rdname color_palettes
#' @export
colorPalette_12 <- c(
  "#9F0162", "#009F81", "#FF5AAF", "#00FCCF", "#8400CD", "#008DF9",
  "#00C2F9", "#FFB2FD", "#A40122", "#E20134", "#FF6E3A", "#FFC33B"
)

#' @rdname color_palettes
#' @export
colorPalette_15 <- c(
  "#68023F", "#008169", "#EF0096", "#00DCB5", "#FFCFE2",
  "#003C86", "#9400E6", "#009FFA", "#FF71FD", "#7CFFFA",
  "#6A0213", "#008607", "#F60239", "#00E307", "#FFDC3D"
)
