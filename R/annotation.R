# Metadata annotation engine -----------------------------------------------------
#
# These internal helpers draw metadata bands onto a rainbow plot at coordinates
# outside the data region (sample metadata to the left/right at negative or > n_target
# x; target metadata above/below at > n_sample or negative y), using
# ggnewscale::new_scale_fill() so each band gets an independent fill scale. Sizes are
# expressed in "cell" units: 1 == the size of one rainbow cell (cells are 1 unit apart
# on both axes), so `width`/`height` of 1 matches a cell, 2 is double, 0.8 is 80%
# (centred on the cell).

# Auto-assign colours per metadata column, cycling the colour-blind-friendly palettes
# by the number of distinct values (interpolating colorPalette_15 beyond 15).
.auto_meta_colors <- function(meta, cols) {
  stats::setNames(lapply(cols, function(cc) {
    x <- meta[[cc]]
    lv <- if (is.factor(x)) {
      levels(droplevels(x))
    } else {
      sort(unique(as.character(x[!is.na(x)])))
    }
    n <- length(lv)
    pal <- if (n == 0) {
      character(0)
    } else if (n <= 8) {
      colorPalette_08[seq_len(n)]
    } else if (n <= 12) {
      colorPalette_12[seq_len(n)]
    } else if (n <= 15) {
      colorPalette_15[seq_len(n)]
    } else {
      grDevices::colorRampPalette(colorPalette_15)(n)
    }
    stats::setNames(pal, lv)
  }), cols)
}

# TRUE for pure white / pure black colours (dropped before building a rank ramp).
.is_white_or_black <- function(cols) {
  rgb <- grDevices::col2rgb(cols)
  colSums(rgb == 255) == 3 | colSums(rgb == 0) == 3
}

# Hue (0-1) of each colour, for ordering a palette before interpolating.
.hue_of <- function(cols) {
  grDevices::rgb2hsv(grDevices::col2rgb(cols))["h", ]
}

# Produce `n` haplotype-rank colours from a supplied palette. If `n` exceeds the palette
# size, drop white/black, sort the rest by hue, and interpolate a ramp to `n` colours.
.expand_rank_palette <- function(palette, n) {
  if (n <= length(palette)) return(palette[seq_len(n)])
  keep <- palette[!.is_white_or_black(palette)]
  if (length(keep) < 2) keep <- palette
  keep <- keep[order(.hue_of(keep))]
  grDevices::colorRampPalette(keep)(n)
}

# Merge user-supplied colours over the auto-assigned ones.
.merge_colors <- function(auto, colors) {
  if (!is.null(colors)) {
    for (nm in names(colors)) auto[[nm]] <- colors[[nm]]
  }
  auto
}

# Resolve a per-band setting (e.g. legend ncol) into a list keyed by column. Accepts
# NULL (all default), a scalar (same for all), a named vector (per column), or an
# unnamed vector applied positionally.
.per_col <- function(x, cols) {
  out <- stats::setNames(vector("list", length(cols)), cols)
  if (is.null(x)) return(out)
  if (!is.null(names(x))) {
    for (nm in intersect(names(x), cols)) out[[nm]] <- x[[nm]]
  } else if (length(x) == 1) {
    for (cc in cols) out[[cc]] <- x[[1]]
  } else {
    for (j in seq_along(cols)) if (j <= length(x)) out[[cols[[j]]]] <- x[[j]]
  }
  out
}

# Build the `guide` argument for a metadata band's fill scale.
.band_guide <- function(legend, ncol, nrow) {
  if (!isTRUE(legend)) return("none")
  if (is.null(ncol) && is.null(nrow)) return("legend")
  ggplot2::guide_legend(ncol = ncol, nrow = nrow)
}

# Extract the legend (guide-box) grob from a built ggplot, for exporting on its own.
.extract_legend_grob <- function(p) {
  gt <- ggplot2::ggplotGrob(p)
  if (!any(grepl("guide-box", gt$layout$name))) {
    stop("The plot has no legend to extract.", call. = FALSE)
  }
  gtable::gtable_filter(gt, "guide-box")
}

# Draw the per-sample metadata sidebar (bands to the left or right of the rainbow).
.add_sample_metadata <- function(p, prepped, sample_meta, sample_key, target_key,
                                 cols, side, width, height, gap, plot_gap, colors,
                                 na_color, legend, legend_ncol, legend_nrow, labels,
                                 target_labels, border) {
  side <- match.arg(side, c("left", "right"))
  if (is.null(cols)) cols <- setdiff(names(sample_meta), sample_key)
  bad <- setdiff(cols, names(sample_meta))
  if (length(bad)) {
    stop("Sample metadata column(s) not found: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  legend_ncol <- .per_col(legend_ncol, cols)
  legend_nrow <- .per_col(legend_nrow, cols)

  samp_levels <- levels(factor(dplyr::pull(prepped, dplyr::all_of(sample_key))))
  tgt_levels  <- levels(factor(dplyr::pull(prepped, dplyr::all_of(target_key))))
  n_t <- length(tgt_levels)

  aligned <- tibble::tibble(!!sample_key := samp_levels) %>%
    dplyr::left_join(sample_meta, by = sample_key)
  aligned[[".ypos"]] <- seq_along(samp_levels)

  pal <- .merge_colors(.auto_meta_colors(aligned, cols), colors)

  step <- width + gap
  mids <- numeric(length(cols))
  for (j in seq_along(cols)) {
    cc <- cols[[j]]
    if (side == "left") {
      inner <- (0.5 - plot_gap) - (j - 1) * step
      xmn <- inner - width; xmx <- inner
    } else {
      inner <- (n_t + 0.5 + plot_gap) + (j - 1) * step
      xmn <- inner; xmx <- inner + width
    }
    mids[[j]] <- (xmn + xmx) / 2
    band <- aligned
    band[[".xmin"]] <- xmn
    band[[".xmax"]] <- xmx
    band[[".ymin"]] <- band[[".ypos"]] - height / 2
    band[[".ymax"]] <- band[[".ypos"]] + height / 2
    p <- p +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_rect(
        data = band,
        mapping = ggplot2::aes(
          xmin = .data[[".xmin"]], xmax = .data[[".xmax"]],
          ymin = .data[[".ymin"]], ymax = .data[[".ymax"]],
          fill = .data[[cc]]
        ),
        color = border
      ) +
      ggplot2::scale_fill_manual(
        name = cc, values = pal[[cc]], na.value = na_color,
        guide = .band_guide(legend, legend_ncol[[cc]], legend_nrow[[cc]])
      )
  }

  brks <- numeric(0); labs <- character(0)
  if (labels)        { brks <- c(brks, mids);          labs <- c(labs, cols) }
  if (target_labels) { brks <- c(brks, seq_len(n_t));  labs <- c(labs, tgt_levels) }
  p +
    ggplot2::scale_x_continuous(
      breaks = if (length(brks)) brks else NULL,
      labels = if (length(brks)) labs else ggplot2::waiver(),
      expand = c(0, 0)
    ) +
    ggplot2::theme(axis.text.x = if (length(brks)) {
      ggplot2::element_text(family = "mono", angle = -90, hjust = 0)
    } else {
      ggplot2::element_blank()
    })
}

# Draw the per-target annotation strip (bands above or below the rainbow).
.add_target_annotation <- function(p, prepped, target_meta, sample_key, target_key,
                                   cols, position, width, height, gap, plot_gap,
                                   colors, na_color, legend, legend_ncol, legend_nrow,
                                   labels, sample_labels, border) {
  position <- match.arg(position, c("top", "bottom"))
  if (is.null(cols)) cols <- setdiff(names(target_meta), target_key)
  bad <- setdiff(cols, names(target_meta))
  if (length(bad)) {
    stop("Target metadata column(s) not found: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  legend_ncol <- .per_col(legend_ncol, cols)
  legend_nrow <- .per_col(legend_nrow, cols)

  samp_levels <- levels(factor(dplyr::pull(prepped, dplyr::all_of(sample_key))))
  tgt_levels  <- levels(factor(dplyr::pull(prepped, dplyr::all_of(target_key))))
  n_s <- length(samp_levels)

  aligned <- tibble::tibble(!!target_key := tgt_levels) %>%
    dplyr::left_join(target_meta, by = target_key)
  aligned[[".xpos"]] <- seq_along(tgt_levels)

  pal <- .merge_colors(.auto_meta_colors(aligned, cols), colors)

  step <- height + gap
  mids <- numeric(length(cols))
  for (j in seq_along(cols)) {
    cc <- cols[[j]]
    if (position == "top") {
      inner <- (n_s + 0.5 + plot_gap) + (j - 1) * step
      ymn <- inner; ymx <- inner + height
    } else {
      inner <- (0.5 - plot_gap) - (j - 1) * step
      ymn <- inner - height; ymx <- inner
    }
    mids[[j]] <- (ymn + ymx) / 2
    band <- aligned
    band[[".ymin"]] <- ymn
    band[[".ymax"]] <- ymx
    band[[".xmin"]] <- band[[".xpos"]] - width / 2
    band[[".xmax"]] <- band[[".xpos"]] + width / 2
    p <- p +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_rect(
        data = band,
        mapping = ggplot2::aes(
          xmin = .data[[".xmin"]], xmax = .data[[".xmax"]],
          ymin = .data[[".ymin"]], ymax = .data[[".ymax"]],
          fill = .data[[cc]]
        ),
        color = border
      ) +
      ggplot2::scale_fill_manual(
        name = cc, values = pal[[cc]], na.value = na_color,
        guide = .band_guide(legend, legend_ncol[[cc]], legend_nrow[[cc]])
      )
  }

  brks <- numeric(0); labs <- character(0)
  if (labels)        { brks <- c(brks, mids);         labs <- c(labs, cols) }
  if (sample_labels) { brks <- c(brks, seq_len(n_s)); labs <- c(labs, samp_levels) }
  p +
    ggplot2::scale_y_continuous(
      breaks = if (length(brks)) brks else NULL,
      labels = if (length(brks)) labs else ggplot2::waiver(),
      expand = c(0, 0)
    ) +
    ggplot2::theme(axis.text.y = if (length(brks)) {
      ggplot2::element_text(family = "mono")
    } else {
      ggplot2::element_blank()
    })
}
