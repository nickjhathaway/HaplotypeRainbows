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

# Reorder colours so that consecutive entries are as visually distinct as possible
# (greedy: each next colour is the farthest, in RGB space, from the previous one). This
# stops a smooth ramp from giving the most abundant ranks (1, 2, 3, ...) near-identical
# colours.
.distinct_order <- function(cols) {
  n <- length(cols)
  if (n <= 2) return(cols)
  d <- as.matrix(stats::dist(t(grDevices::col2rgb(cols))))
  ord <- integer(n)
  used <- logical(n)
  ord[1] <- 1L
  used[1] <- TRUE
  for (i in 2:n) {
    cand <- which(!used)
    ord[i] <- cand[which.max(d[ord[i - 1], cand])]
    used[ord[i]] <- TRUE
  }
  cols[ord]
}

# Produce `n` haplotype-rank colours from a supplied palette. If `n` exceeds the palette
# size, drop white/black, sort the rest by hue, interpolate a ramp to `n` colours, then
# reorder so adjacent ranks are maximally distinct.
.expand_rank_palette <- function(palette, n) {
  if (n <= length(palette)) return(palette[seq_len(n)])
  keep <- palette[!.is_white_or_black(palette)]
  if (length(keep) < 2) keep <- palette
  keep <- keep[order(.hue_of(keep))]
  .distinct_order(grDevices::colorRampPalette(keep)(n))
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

# Coerce the given metadata columns to factors, honouring any user-supplied level
# order (per column, via a named list). Supplied levels come first; any remaining
# values are appended in their natural (sorted) order. This controls both the colour
# assignment and the legend entry order.
.apply_level_order <- function(aligned, cols, level_order) {
  lo <- .per_col(level_order, cols)
  for (cc in cols) {
    x <- as.character(aligned[[cc]])
    present <- sort(unique(x[!is.na(x)]))
    ord <- lo[[cc]]
    full <- if (is.null(ord)) {
      present
    } else {
      c(as.character(ord), setdiff(present, as.character(ord)))
    }
    aligned[[cc]] <- factor(x, levels = full)
  }
  aligned
}

# Levels of a (possibly spacer-augmented) factor column, WITHOUT dropping unused
# levels — so spacer levels inserted between clusters keep their y positions.
.factor_levels <- function(prepped, key) {
  x <- dplyr::pull(prepped, dplyr::all_of(key))
  if (is.factor(x)) levels(x) else levels(factor(x))
}

# Build the `guide` argument for a metadata band's fill scale.
.band_guide <- function(legend, ncol, nrow) {
  if (!isTRUE(legend)) return("none")
  if (is.null(ncol) && is.null(nrow)) return("legend")
  ggplot2::guide_legend(ncol = ncol, nrow = nrow)
}

# Descend a guide-box to the nested gtable that actually holds the legend content
# (the one with real absolute width). The outer guide-box has flexible `null` spacer
# cells that would expand — and mis-place the legend — if drawn on its own.
.legend_content_grob <- function(guide_box) {
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf, width = 50, height = 50)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)
  abs_width <- function(g) {
    w <- tryCatch(grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
                  error = function(e) 0)
    if (!is.finite(w)) 0 else w
  }
  best <- guide_box
  best_w <- abs_width(guide_box)
  recurse <- function(g) {
    if (!is.null(g$grobs)) {
      for (ch in g$grobs) {
        if (inherits(ch, "gtable")) {
          cw <- abs_width(ch)
          if (cw > best_w) {
            best <<- ch
            best_w <<- cw
          }
          recurse(ch)
        }
      }
    }
  }
  recurse(guide_box)
  best
}

# Extract the legend grob from a built ggplot, for exporting/composing on its own.
.extract_legend_grob <- function(p) {
  gt <- ggplot2::ggplotGrob(p)
  if (!any(grepl("guide-box", gt$layout$name))) {
    stop("The plot has no legend to extract.", call. = FALSE)
  }
  .legend_content_grob(gtable::gtable_filter(gt, "guide-box"))
}

# Estimate the natural size (inches) of a legend grob, so it can be exported without
# clipping. The guide-box's own widths can be a flexible (null) cell whose real content
# lives in a nested gtable, so we recurse and take the largest absolute size found.
# Measured on a throwaway device so text-based (grobwidth) units resolve.
.legend_grob_size <- function(legend, margin = 0.2) {
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf, width = 50, height = 50)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)
  measure <- function(g) {
    w <- tryCatch(grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
                  error = function(e) 0)
    h <- tryCatch(grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE),
                  error = function(e) 0)
    if (!is.finite(w)) w <- 0
    if (!is.finite(h)) h <- 0
    if (!is.null(g$grobs)) {
      for (ch in g$grobs) {
        if (inherits(ch, "gtable")) {
          s <- measure(ch)
          w <- max(w, s[[1]])
          h <- max(h, s[[2]])
        }
      }
    }
    c(w, h)
  }
  sz <- measure(legend)
  w <- if (sz[[1]] > 0) sz[[1]] else 6
  h <- if (sz[[2]] > 0) sz[[2]] else 8
  c(width = w + 2 * margin, height = h + 2 * margin)
}

# Draw alternating translucent background bands, one per cluster group, spanning the
# plot width. `groups` is a tibble with the sample column and a `cluster` factor.
# `expand` pushes the band edges `expand` cells past the plotting cells on both sides
# (widening the panel), so the bands show as clearly visible strips in the margin even
# when the cell grid is dense.
.add_cluster_bands <- function(p, prepped, sample_key, target_key, groups, colors,
                               extend_left, extend_right, expand, border) {
  samp_levels <- .factor_levels(prepped, sample_key)
  tgt_levels  <- .factor_levels(prepped, target_key)
  n_t <- length(tgt_levels)
  pos <- stats::setNames(seq_along(samp_levels), samp_levels)

  g <- as.data.frame(groups)
  g[["pos"]] <- pos[as.character(g[[sample_key]])]
  g <- g[!is.na(g[["pos"]]), , drop = FALSE]
  bands <- g %>%
    dplyr::group_by(.data[["cluster"]]) %>%
    dplyr::summarise(ymin = min(.data[["pos"]]) - 0.5,
                     ymax = max(.data[["pos"]]) + 0.5, .groups = "drop") %>%
    dplyr::arrange(.data[["ymin"]]) %>%
    dplyr::mutate(fill = colors[(dplyr::row_number() - 1) %% length(colors) + 1])

  xmn <- 0.5 - extend_left - expand
  xmx <- n_t + 0.5 + extend_right + expand
  p +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(
      data = bands,
      mapping = ggplot2::aes(
        xmin = xmn, xmax = xmx,
        ymin = .data[["ymin"]], ymax = .data[["ymax"]], fill = .data[["fill"]]
      ),
      colour = border
    ) +
    ggplot2::scale_fill_identity()
}

# Draw the per-sample metadata sidebar (bands to the left or right of the rainbow).
.add_sample_metadata <- function(p, prepped, sample_meta, sample_key, target_key,
                                 cols, side, width, height, gap, plot_gap, colors,
                                 na_color, legend, legend_ncol, legend_nrow,
                                 level_order, labels, target_labels, border) {
  side <- match.arg(side, c("left", "right"))
  if (is.null(cols)) cols <- setdiff(names(sample_meta), sample_key)
  bad <- setdiff(cols, names(sample_meta))
  if (length(bad)) {
    stop("Sample metadata column(s) not found: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  legend_ncol <- .per_col(legend_ncol, cols)
  legend_nrow <- .per_col(legend_nrow, cols)

  samp_levels <- .factor_levels(prepped, sample_key)
  tgt_levels  <- .factor_levels(prepped, target_key)
  n_t <- length(tgt_levels)

  aligned <- tibble::tibble(!!sample_key := samp_levels) %>%
    dplyr::left_join(sample_meta, by = sample_key)
  aligned[[".ypos"]] <- seq_along(samp_levels)
  # drop spacer rows (levels with no data, e.g. inter-cluster gaps) so we don't draw
  # empty metadata boxes there; their positions are already baked into .ypos.
  real <- unique(as.character(dplyr::pull(prepped, dplyr::all_of(sample_key))))
  aligned <- aligned[aligned[[sample_key]] %in% real, , drop = FALSE]
  aligned <- .apply_level_order(aligned, cols, level_order)

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
                                   level_order, labels, sample_labels, border) {
  position <- match.arg(position, c("top", "bottom"))
  if (is.null(cols)) cols <- setdiff(names(target_meta), target_key)
  bad <- setdiff(cols, names(target_meta))
  if (length(bad)) {
    stop("Target metadata column(s) not found: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  legend_ncol <- .per_col(legend_ncol, cols)
  legend_nrow <- .per_col(legend_nrow, cols)

  samp_levels <- .factor_levels(prepped, sample_key)
  tgt_levels  <- .factor_levels(prepped, target_key)
  n_s <- length(samp_levels)

  aligned <- tibble::tibble(!!target_key := tgt_levels) %>%
    dplyr::left_join(target_meta, by = target_key)
  aligned[[".xpos"]] <- seq_along(tgt_levels)
  aligned <- .apply_level_order(aligned, cols, level_order)

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
