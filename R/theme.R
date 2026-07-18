# Internal plotting helpers ----------------------------------------------------

# The base theme shared by all rainbow plots (formerly the inline "sofonias_theme").
.rainbow_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.3),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.3),
      text = ggplot2::element_text(size = 12, family = "Helvetica"),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(size = 12),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

# Build the y-axis breaks/labels for the sample factor. Samples are drawn at their
# integer factor position; gaps (missing integer levels) get blank labels.
.sample_axis <- function(prep_data, sample_col) {
  sample_names <- prep_data %>%
    dplyr::group_by() %>%
    dplyr::select(dplyr::all_of(sample_col)) %>%
    unique() %>%
    dplyr::arrange(.data[[sample_col]])
  vals <- sample_names[[1]]
  min_level <- min(as.numeric(vals))
  max_level <- max(as.numeric(vals))
  labels_by_level <- as.character(vals)
  names(labels_by_level) <- as.numeric(vals)
  plot_labels <- vapply(min_level:max_level, function(lev) {
    key <- as.character(lev)
    if (key %in% names(labels_by_level)) labels_by_level[[key]] else ""
  }, character(1))
  list(breaks = min_level:max_level, labels = plot_labels)
}
