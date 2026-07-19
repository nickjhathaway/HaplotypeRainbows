# Internal plotting helpers ----------------------------------------------------

# The base theme shared by all rainbow plots (formerly the inline "sofonias_theme").
.rainbow_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.3),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.3),
      text = ggplot2::element_text(size = 12, family = "Helvetica"),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(size = 12),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}
