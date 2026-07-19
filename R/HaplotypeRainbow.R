#' HaplotypeRainbow
#'
#' An R6 class for building haplotype "rainbow" plots. The class carries the mapping
#' of your four data columns (sample, target, haplotype identifier and relative counts)
#' so it only has to be supplied once. Methods that transform state (`prep()`,
#' `sort_by_clustering()`, `set_sample_order()`, `sort_alphabetical()`) mutate the object
#' and return it invisibly, so they can be chained; `plot()` returns a [ggplot2::ggplot]
#' object.
#'
#' Column defaults follow the Portable Microhaplotype Object (PMO) convention:
#' `library_sample_name`, `target_name`, `seq` and `reads`.
#'
#' @examples
#' \dontrun{
#' rb <- HaplotypeRainbow$new(my_allele_table)
#' rb$prep(sort = "population")$sort_by_clustering()
#' rb$plot()
#' }
#'
#' @export
HaplotypeRainbow <- R6::R6Class(
  "HaplotypeRainbow",
  public = list(

    #' @description Create a new HaplotypeRainbow.
    #' @param data A data frame / tibble with one row per (sample, target, haplotype).
    #' @param sample_col Name of the sample identifier column.
    #' @param target_col Name of the target / locus column.
    #' @param popuid_col Name of the within-target haplotype identifier column.
    #' @param rel_abund_col Name of the relative counts column (raw counts are fine;
    #'   they are normalised to within-sample/target fractions internally).
    initialize = function(data,
                          sample_col    = "library_sample_name",
                          target_col    = "target_name",
                          popuid_col    = "seq",
                          rel_abund_col = "reads") {
      stopifnot(is.data.frame(data))
      cols <- list(sample = sample_col, target = target_col,
                   popuid = popuid_col, rel_abund = rel_abund_col)
      missing <- setdiff(unlist(cols), colnames(data))
      if (length(missing)) {
        stop("Column(s) not found in data: ", paste(missing, collapse = ", "),
             call. = FALSE)
      }
      private$data <- data
      private$cols <- cols
      invisible(self)
    },

    #' @description Prep the data for a rainbow plot (haplotype colours rotate across
    #'   targets). Stores the prepped table internally.
    #' @param sort Haplotype ordering within each cell: "population_rank" (order by
    #'   population rank) or "within_sample_freq" (order by within-sample fraction).
    #' @param min_pop_size Drop targets with fewer than this many unique haplotypes.
    #' @param color_period Number of hue steps to rotate across targets.
    #' @param bar_height Height of the full stacked bar per sample; < 1 leaves a gap.
    #' @param mark_invariant Flag single-haplotype (invariant) targets so that
    #'   `plot(rank_colors = TRUE)` can colour them separately (black by default).
    #' @return The object, invisibly (chainable).
    prep = function(sort = c("population_rank", "within_sample_freq"),
                    min_pop_size = 1, color_period = 11, bar_height = 0.80,
                    mark_invariant = FALSE) {
      sort <- match.arg(sort)
      private$prepped <- .prep_engine(
        private$data, private$cols, sort = sort,
        min_pop_size = min_pop_size, color_period = color_period,
        bar_height = bar_height
      )
      if (isTRUE(mark_invariant)) {
        private$prepped <- private$prepped %>%
          dplyr::mutate(is_invariant = maxPopid == 1)
      }
      private$sort_mode <- sort
      private$color_period <- color_period
      invisible(self)
    },

    #' @description Prep the data for a shade plot (haplotypes shaded from a per-target
    #'   base colour instead of a rotating rainbow). Stores the prepped table internally.
    #' @param min_pop_size Drop targets with fewer than this many unique haplotypes.
    #' @param bar_height Height of the full stacked bar per sample; < 1 leaves a gap.
    #' @param base_colors Base colours to generate the per-target shades from.
    #' @return The object, invisibly (chainable).
    prep_shade = function(min_pop_size = 1, bar_height = 0.80,
                          base_colors = c("#e41a1c", "#377eb8", "#4daf4a",
                                          "#984ea3", "#ff7f00", "#ffff33")) {
      private$prepped <- .prep_engine(
        private$data, private$cols, sort = "shade",
        min_pop_size = min_pop_size, bar_height = bar_height,
        base_colors = base_colors
      )
      private$sort_mode <- "shade"
      invisible(self)
    },

    #' @description Reorder samples so that similar samples (by haplotype sharing) sit
    #'   next to each other, via hierarchical clustering (ward.D2).
    #' @param by_major_allele Cluster on only the major haplotype per (sample, target).
    #' @param coverage_cutoff Keep only targets present in at least this fraction of
    #'   samples for the clustering.
    #' @return The object, invisibly (chainable).
    sort_by_clustering = function(by_major_allele = FALSE, coverage_cutoff = 0.80) {
      private$require_prepped()
      private$prepped <- private$resort_clustering(
        private$prepped, private$cols, coverage_cutoff, by_major_allele
      )
      invisible(self)
    },

    #' @description Set the sample order explicitly.
    #' @param levels Character vector of sample names in the desired order.
    #' @return The object, invisibly (chainable).
    set_sample_order = function(levels) {
      private$require_prepped()
      sc <- private$cols$sample
      present <- as.character(unique(dplyr::pull(private$prepped, dplyr::all_of(sc))))
      missing <- setdiff(present, as.character(levels))
      ordered <- c(as.character(levels), missing)
      private$prepped <- private$prepped %>%
        dplyr::mutate("{sc}" := factor(.data[[sc]], levels = ordered))
      invisible(self)
    },

    #' @description Order samples alphabetically (the default).
    #' @return The object, invisibly (chainable).
    sort_alphabetical = function() {
      private$require_prepped()
      sc <- private$cols$sample
      lv <- sort(as.character(unique(dplyr::pull(private$prepped, dplyr::all_of(sc)))))
      private$prepped <- private$prepped %>%
        dplyr::mutate("{sc}" := factor(.data[[sc]], levels = lv))
      invisible(self)
    },

    #' @description Attach per-sample metadata to the object (used by the metadata
    #'   sidebar and by `sort_samples_by_meta()`).
    #' @param meta A data frame of sample metadata.
    #' @param match_col Name of the column in `meta` that holds the sample identifiers.
    #' @param cols Metadata columns to keep (default: all columns except `match_col`).
    #' @param add If `FALSE` (default) replace any existing sample metadata; if `TRUE`
    #'   merge these columns onto the existing sample metadata (joined on the samples).
    #' @return The object, invisibly (chainable).
    set_sample_meta = function(meta, match_col, cols = NULL, add = FALSE) {
      private$set_meta("sample", meta, match_col, cols, add)
    },

    #' @description Attach per-target metadata to the object (used by the target
    #'   annotation strip and by `sort_targets_by_meta()`).
    #' @param meta A data frame of target metadata.
    #' @param match_col Name of the column in `meta` that holds the target names.
    #' @param cols Metadata columns to keep (default: all columns except `match_col`).
    #' @param add If `FALSE` (default) replace any existing target metadata; if `TRUE`
    #'   merge these columns onto the existing target metadata (joined on the targets).
    #' @return The object, invisibly (chainable).
    set_target_meta = function(meta, match_col, cols = NULL, add = FALSE) {
      private$set_meta("target", meta, match_col, cols, add)
    },

    #' @description Return the stored sample metadata (or `NULL`).
    #' @return A data frame.
    get_sample_meta = function() private$sample_meta,

    #' @description Return the stored target metadata (or `NULL`).
    #' @return A data frame.
    get_target_meta = function() private$target_meta,

    #' @description Order samples by one or more sample-metadata columns. Samples with
    #'   no metadata are placed at the end.
    #' @param cols Character vector of sample-metadata column names to order by.
    #' @param desc Sort descending instead of ascending.
    #' @return The object, invisibly (chainable).
    sort_samples_by_meta = function(cols, desc = FALSE) {
      private$require_prepped()
      private$prepped <- private$sort_by_meta("sample", cols, desc)
      invisible(self)
    },

    #' @description Order targets by one or more target-metadata columns. Targets with
    #'   no metadata are placed at the end.
    #' @param cols Character vector of target-metadata column names to order by.
    #' @param desc Sort descending instead of ascending.
    #' @return The object, invisibly (chainable).
    sort_targets_by_meta = function(cols, desc = FALSE) {
      private$require_prepped()
      private$prepped <- private$sort_by_meta("target", cols, desc)
      invisible(self)
    },

    #' @description Set the target (column) order explicitly.
    #' @param levels Character vector of target names in the desired order.
    #' @return The object, invisibly (chainable).
    set_target_order = function(levels) {
      private$require_prepped()
      tc <- private$cols$target
      present <- levels(factor(dplyr::pull(private$prepped, dplyr::all_of(tc))))
      ordered <- c(as.character(levels),
                   setdiff(present, as.character(levels)))
      private$prepped <- private$prepped %>%
        dplyr::mutate("{tc}" := factor(.data[[tc]], levels = ordered))
      invisible(self)
    },

    #' @description Build the rainbow plot.
    #' @param style "rainbow" (gradient over rotating hues) or "shade" (identity colours
    #'   from shade prep). Defaults to "shade" when the data was prepped with
    #'   `sort = "shade"`, otherwise "rainbow".
    #' @param colors Colour ramp for rainbow style. Length should match `color_period`.
    #' @param color_col Column to map to fill for rainbow style
    #'   ("popidFracLogColor" or "popidFracRegColor").
    #' @param shade_col Identity-colour column for shade style.
    #' @param rank_colors Colour haplotypes by a discrete rank (with a legend) instead
    #'   of the continuous rainbow gradient. When there are more ranks than palette
    #'   colours, a colour ramp is interpolated from the palette (with a warning);
    #'   invariant targets (see `prep(mark_invariant = TRUE)`) get `invariant_color`.
    #'   Replaces the old harvest-fills-from-`ggplot_build` recipe.
    #' @param rank_palette Colours for haplotype ranks 1, 2, ... (rank colouring); if
    #'   there are more ranks than colours, a ramp is interpolated from these.
    #' @param invariant_color Colour for invariant (single-haplotype) targets.
    #' @param rank_legend_title Legend title used with `rank_colors`.
    #' @param x_axis_labels Show target names on the x-axis.
    #' @param y_axis_labels Show sample names on the y-axis.
    #' @return A [ggplot2::ggplot] object.
    plot = function(style = NULL,
                    colors = RColorBrewer::brewer.pal(11, "Spectral"),
                    color_col = "popidFracLogColor",
                    shade_col = "h_color_byFreq_mod",
                    rank_colors = FALSE,
                    rank_palette = colorPalette_12,
                    invariant_color = "#000000",
                    rank_legend_title = "Microhaplotype Rank",
                    x_axis_labels = TRUE, y_axis_labels = TRUE) {
      private$require_prepped()
      if (is.null(style)) {
        style <- if (identical(private$sort_mode, "shade")) "shade" else "rainbow"
      }
      style <- match.arg(style, c("rainbow", "shade"))
      if (style == "rainbow" && identical(private$sort_mode, "shade")) {
        stop("style = 'rainbow' needs data prepped with $prep().", call. = FALSE)
      }
      if (style == "shade" && !identical(private$sort_mode, "shade")) {
        stop("style = 'shade' needs data prepped with $prep_shade().", call. = FALSE)
      }
      if (isTRUE(rank_colors) && identical(private$sort_mode, "shade")) {
        stop("rank_colors requires data prepped with $prep() (not $prep_shade()).",
             call. = FALSE)
      }
      rank <- if (isTRUE(rank_colors)) {
        list(palette = rank_palette, invariant_color = invariant_color,
             title = rank_legend_title)
      } else {
        NULL
      }
      fill_col <- if (style == "rainbow") color_col else shade_col
      private$build_plot(style, fill_col, colors, x_axis_labels, y_axis_labels, rank)
    },

    #' @description Return the discrete rank -> colour map used by
    #'   `plot(rank_colors = TRUE)` (haplotype ranks plus the invariant colour). When
    #'   there are more ranks than palette colours a ramp is interpolated.
    #' @param n_ranks Number of ranks to produce colours for (default: the maximum
    #'   haplotype rank in the prepped data, or the palette length if not prepped).
    #' @param rank_palette Colours for ranks 1, 2, ...
    #' @param invariant_color Colour for invariant targets.
    #' @return A named character vector.
    get_rank_colors = function(n_ranks = NULL, rank_palette = colorPalette_12,
                               invariant_color = "#000000") {
      if (is.null(n_ranks)) {
        n_ranks <- if (!is.null(private$prepped) &&
                       "popid" %in% names(private$prepped)) {
          max(private$prepped[["popid"]], na.rm = TRUE)
        } else {
          length(rank_palette)
        }
      }
      pal <- .expand_rank_palette(rank_palette, n_ranks)
      c(stats::setNames(pal, as.character(seq_len(n_ranks))),
        invariant = invariant_color)
    },

    #' @description Add a per-sample metadata sidebar (coloured bands) to a rainbow
    #'   plot, using the metadata attached with `set_sample_meta()`.
    #' @param p A ggplot returned by `plot()`.
    #' @param cols Metadata columns to draw (default: all in the sample-metadata slot),
    #'   one band each, in order outward from the plot.
    #' @param side "left" or "right" of the rainbow.
    #' @param width Band width in cell units (1 = one rainbow cell wide).
    #' @param height Band height as a fraction of a cell (1 = full cell, centred).
    #' @param gap Gap between bands in cell units (0 = flush).
    #' @param plot_gap Gap in cell units between the rainbow cells and the first band
    #'   (0 = flush against the plot).
    #' @param colors Optional named list (per column) of value -> colour to override
    #'   the auto-assigned colours.
    #' @param na_color Fill for samples missing that metadata value.
    #' @param legend Draw a colour legend for each band.
    #' @param legend_ncol,legend_nrow Number of columns / rows for the band legends,
    #'   to keep tall legends on the page. A scalar applies to every band; a named
    #'   vector (e.g. `c(country = 3)`) sets it per column.
    #' @param labels Draw the band (column-name) labels on the x-axis. Independent of
    #'   `target_labels`, so band labels can stay on when target names are off.
    #' @param target_labels Keep the target names on the x-axis.
    #' @param border Border colour for the band rectangles.
    #' @return A [ggplot2::ggplot] object.
    add_sample_metadata = function(p, cols = NULL, side = "left", width = 1,
                                   height = 1, gap = 0, plot_gap = 0, colors = NULL,
                                   na_color = "grey80", legend = TRUE,
                                   legend_ncol = NULL, legend_nrow = NULL,
                                   labels = TRUE, target_labels = TRUE,
                                   border = "black") {
      private$require_prepped()
      if (is.null(private$sample_meta)) {
        stop("No sample metadata set; call $set_sample_meta() first.", call. = FALSE)
      }
      .add_sample_metadata(p, private$prepped, private$sample_meta,
                           private$cols$sample, private$cols$target, cols, side,
                           width, height, gap, plot_gap, colors, na_color, legend,
                           legend_ncol, legend_nrow, labels, target_labels, border)
    },

    #' @description Add a per-target annotation strip (coloured bands) above or below a
    #'   rainbow plot, using the metadata attached with `set_target_meta()`.
    #' @param p A ggplot returned by `plot()`.
    #' @param cols Metadata columns to draw (default: all in the target-metadata slot),
    #'   one band each, stacked outward from the plot.
    #' @param position "top" or "bottom" of the rainbow.
    #' @param width Band width as a fraction of a cell (1 = full cell, centred).
    #' @param height Band thickness in cell units (1 = one rainbow cell tall).
    #' @param gap Gap between bands in cell units (0 = flush).
    #' @param plot_gap Gap in cell units between the rainbow cells and the first band
    #'   (0 = flush against the plot).
    #' @param colors Optional named list (per column) of value -> colour to override
    #'   the auto-assigned colours.
    #' @param na_color Fill for targets missing that metadata value.
    #' @param legend Draw a colour legend for each band.
    #' @param legend_ncol,legend_nrow Number of columns / rows for the band legends.
    #'   A scalar applies to every band; a named vector sets it per column.
    #' @param labels Draw the band (column-name) labels on the y-axis. Independent of
    #'   `sample_labels`.
    #' @param sample_labels Keep the sample names on the y-axis.
    #' @param border Border colour for the band rectangles.
    #' @return A [ggplot2::ggplot] object.
    add_target_annotation = function(p, cols = NULL, position = "top", width = 1,
                                     height = 1, gap = 0, plot_gap = 0, colors = NULL,
                                     na_color = "grey80", legend = TRUE,
                                     legend_ncol = NULL, legend_nrow = NULL,
                                     labels = TRUE, sample_labels = TRUE,
                                     border = "black") {
      private$require_prepped()
      if (is.null(private$target_meta)) {
        stop("No target metadata set; call $set_target_meta() first.", call. = FALSE)
      }
      .add_target_annotation(p, private$prepped, private$target_meta,
                             private$cols$sample, private$cols$target, cols, position,
                             width, height, gap, plot_gap, colors, na_color, legend,
                             legend_ncol, legend_nrow, labels, sample_labels, border)
    },

    #' @description Suggested figure dimensions (inches) for the current data, scaling
    #'   width with the number of targets and height with the number of samples.
    #' @param cell_width Inches per target column.
    #' @param cell_height Inches per sample row.
    #' @param min_width Minimum width.
    #' @param min_height Minimum height.
    #' @param extra_width Padding added to width (e.g. for a sidebar / legends).
    #' @param extra_height Padding added to height (e.g. for a target strip / legends).
    #' @return A list with `width` and `height` (inches).
    dims = function(cell_width = 0.3, cell_height = 0.3, min_width = 6,
                    min_height = 6, extra_width = 0, extra_height = 0) {
      private$require_prepped()
      n_t <- dplyr::n_distinct(
        dplyr::pull(private$prepped, dplyr::all_of(private$cols$target))
      )
      n_s <- dplyr::n_distinct(
        dplyr::pull(private$prepped, dplyr::all_of(private$cols$sample))
      )
      list(width  = max(min_width,  n_t * cell_width)  + extra_width,
           height = max(min_height, n_s * cell_height) + extra_height)
    },

    #' @description Save a plot to PDF, sized automatically from the data when `width` /
    #'   `height` are not supplied. Uses `cairo_pdf` by default (which keeps hyphens as
    #'   hyphens rather than converting them to en/em dashes); pass `device = "pdf"` for
    #'   the base device (e.g. on Windows, where the cairo device can misbehave).
    #' @param p A ggplot (e.g. from `plot()` or `add_sample_metadata()`).
    #' @param file Output file path.
    #' @param width Width in inches (default from `dims()`).
    #' @param height Height in inches (default from `dims()`).
    #' @param device "cairo" (default) or "pdf".
    #' @param ... Passed to the graphics device.
    #' @return The file path, invisibly.
    save_pdf = function(p, file, width = NULL, height = NULL,
                        device = c("cairo", "pdf"), ...) {
      device <- match.arg(device)
      if (is.null(width) || is.null(height)) {
        d <- self$dims()
        if (is.null(width)) width <- d$width
        if (is.null(height)) height <- d$height
      }
      if (device == "cairo") {
        grDevices::cairo_pdf(file, width = width, height = height, ...)
      } else {
        grDevices::pdf(file, width = width, height = height,
                       useDingbats = FALSE, ...)
      }
      on.exit(grDevices::dev.off())
      print(p)
      invisible(file)
    },

    #' @description Return a copy of a plot with all legends removed. Combine with
    #'   `save_legend_pdf()` to export the plot and its legend separately (e.g. to
    #'   place the legend elsewhere in a figure).
    #' @param p A ggplot.
    #' @return A [ggplot2::ggplot] object with `legend.position = "none"`.
    drop_legends = function(p) {
      p + ggplot2::theme(legend.position = "none")
    },

    #' @description Extract the legend (guide-box) of a plot as a grob, for composing
    #'   it elsewhere (e.g. with patchwork / cowplot / grid).
    #' @param p A ggplot.
    #' @return A grid grob (gtable).
    extract_legend = function(p) {
      .extract_legend_grob(p)
    },

    #' @description Save just the legend of a plot to PDF (so it can be combined in
    #'   post with the plot exported via `save_pdf(drop_legends(p), ...)`).
    #' @param p A ggplot.
    #' @param file Output file path.
    #' @param width Width in inches.
    #' @param height Height in inches.
    #' @param device "cairo" (default) or "pdf".
    #' @param ... Passed to the graphics device.
    #' @return The file path, invisibly.
    save_legend_pdf = function(p, file, width = 6, height = 8,
                               device = c("cairo", "pdf"), ...) {
      device <- match.arg(device)
      legend <- .extract_legend_grob(p)
      if (device == "cairo") {
        grDevices::cairo_pdf(file, width = width, height = height, ...)
      } else {
        grDevices::pdf(file, width = width, height = height,
                       useDingbats = FALSE, ...)
      }
      on.exit(grDevices::dev.off())
      grid::grid.newpage()
      grid::grid.draw(legend)
      invisible(file)
    },

    #' @description Return the prepped data frame (or `NULL` if `prep()` not yet called).
    #' @return A data frame.
    get_prepped = function() private$prepped,

    #' @description Print a short summary.
    #' @param ... Ignored.
    print = function(...) {
      cat("<HaplotypeRainbow>\n")
      cat("  columns: sample =", private$cols$sample,
          "| target =", private$cols$target,
          "| haplotype =", private$cols$popuid,
          "| counts =", private$cols$rel_abund, "\n")
      cat("  rows in:", nrow(private$data), "\n")
      if (is.null(private$prepped)) {
        cat("  prepped: <not yet - call $prep()>\n")
      } else {
        cat("  prepped:", nrow(private$prepped), "rows | sort =",
            private$sort_mode, "\n")
      }
      invisible(self)
    }
  ),

  private = list(
    data = NULL,
    cols = NULL,
    prepped = NULL,
    sort_mode = NULL,
    color_period = 11,
    sample_meta = NULL,
    target_meta = NULL,

    require_prepped = function() {
      if (is.null(private$prepped)) {
        stop("Call $prep() before this operation.", call. = FALSE)
      }
    },

    # Store sample/target metadata on the object. `kind` is "sample" or "target";
    # the meta's `match_col` is renamed to the data's key column so downstream joins
    # are uniform. Errors on missing columns; warns + NA-fills entities absent from meta.
    set_meta = function(kind, meta, match_col, cols, add) {
      stopifnot(is.data.frame(meta))
      key_out <- if (kind == "sample") private$cols$sample else private$cols$target
      if (!match_col %in% names(meta)) {
        stop("match_col '", match_col, "' not found in the metadata table.",
             call. = FALSE)
      }
      if (is.null(cols)) cols <- setdiff(names(meta), match_col)
      missing_cols <- setdiff(cols, names(meta))
      if (length(missing_cols)) {
        stop("Metadata column(s) not found: ",
             paste(missing_cols, collapse = ", "), call. = FALSE)
      }
      new_meta <- as.data.frame(meta)[, c(match_col, cols), drop = FALSE]
      names(new_meta)[1] <- key_out
      new_meta[[key_out]] <- as.character(new_meta[[key_out]])
      new_meta <- new_meta[!duplicated(new_meta[[key_out]]), , drop = FALSE]

      data_ids <- unique(as.character(
        dplyr::pull(private$data, dplyr::all_of(key_out))
      ))
      missing_ids <- setdiff(data_ids, new_meta[[key_out]])
      if (length(missing_ids)) {
        show <- missing_ids[seq_len(min(5L, length(missing_ids)))]
        warning(length(missing_ids), " ", kind,
                "(s) in the data are missing from the metadata (filled with NA): ",
                paste(show, collapse = ", "),
                if (length(missing_ids) > 5L) ", ..." else "", call. = FALSE)
      }

      slot <- paste0(kind, "_meta")
      if (isTRUE(add) && !is.null(private[[slot]])) {
        private[[slot]] <- dplyr::full_join(private[[slot]], new_meta, by = key_out)
      } else {
        private[[slot]] <- new_meta
      }
      invisible(self)
    },

    # Reorder the sample (or target) factor of the prepped data by metadata columns.
    sort_by_meta = function(kind, cols, desc) {
      slot <- paste0(kind, "_meta")
      meta <- private[[slot]]
      if (is.null(meta)) {
        stop("No ", kind, " metadata set; call $set_", kind, "_meta() first.",
             call. = FALSE)
      }
      key <- if (kind == "sample") private$cols$sample else private$cols$target
      bad <- setdiff(cols, setdiff(names(meta), key))
      if (length(bad)) {
        stop(kind, " metadata column(s) not found: ",
             paste(bad, collapse = ", "), call. = FALSE)
      }
      present <- levels(factor(dplyr::pull(private$prepped, dplyr::all_of(key))))
      keys <- lapply(cols, function(cc) {
        s <- rlang::sym(cc)
        if (isTRUE(desc)) rlang::expr(dplyr::desc(!!s)) else s
      })
      ordered <- meta %>%
        dplyr::filter(.data[[key]] %in% present) %>%
        dplyr::arrange(!!!keys) %>%
        dplyr::pull(dplyr::all_of(key))
      ordered <- as.character(ordered)
      ordered <- c(ordered, setdiff(present, ordered))
      private$prepped %>%
        dplyr::mutate("{key}" := factor(.data[[key]], levels = ordered))
    },

    # Port of the historical resort_prepped_samples_by_clustering().
    resort_clustering = function(prepped, cols, coverage_cutoff, by_major_allele) {
      s_sym <- rlang::sym(cols$sample)
      t_sym <- rlang::sym(cols$target)
      h_sym <- rlang::sym(cols$popuid)
      a_sym <- rlang::sym(cols$rel_abund)

      n_samples_total <- dplyr::n_distinct(dplyr::pull(prepped, !!s_sym))

      targets_keep <- prepped %>%
        dplyr::ungroup() %>%
        dplyr::group_by(!!t_sym) %>%
        dplyr::summarise(sample_count = dplyr::n_distinct(!!s_sym), .groups = "drop") %>%
        dplyr::mutate(sample_freq = sample_count / n_samples_total) %>%
        dplyr::filter(sample_freq >= coverage_cutoff) %>%
        dplyr::pull(!!t_sym) %>%
        unique()
      if (length(targets_keep) == 0) {
        stop("No targets above the sample-coverage cutoff.", call. = FALSE)
      }

      dat <- prepped %>% dplyr::filter(!!t_sym %in% targets_keep)
      if (by_major_allele) {
        dat <- dat %>%
          dplyr::group_by(!!s_sym, !!t_sym) %>%
          dplyr::slice_max(order_by = !!a_sym, n = 1, with_ties = FALSE) %>%
          dplyr::ungroup()
      }

      dat_sp <- dat %>%
        dplyr::mutate(new_identifier = paste0(!!t_sym, "-", !!h_sym)) %>%
        dplyr::ungroup() %>%
        dplyr::select(!!s_sym, new_identifier) %>%
        unique() %>%
        dplyr::mutate(marker = 1L) %>%
        tidyr::pivot_wider(names_from = new_identifier, values_from = marker,
                           values_fill = 0L)

      dat_sp_mat <- as.matrix(dat_sp[, 2:ncol(dat_sp)])
      rownames(dat_sp_mat) <- dplyr::pull(dat_sp, !!s_sym)
      hc <- stats::hclust(stats::dist(dat_sp_mat), method = "ward.D2")
      sample_levels <- rownames(dat_sp_mat)[hc$order]

      missing_samples <- prepped %>%
        dplyr::filter(!(!!s_sym %in% sample_levels)) %>%
        dplyr::ungroup() %>%
        dplyr::pull(!!s_sym) %>%
        as.character() %>%
        unique()

      sample_levels <- c(as.character(sample_levels), missing_samples)
      prepped %>%
        dplyr::mutate("{cols$sample}" := factor(!!s_sym, levels = sample_levels))
    },

    build_plot = function(style, fill_col, colors, x_axis_labels, y_axis_labels,
                          rank = NULL) {
      prep_data <- private$prepped
      sc <- private$cols$sample
      tc <- private$cols$target
      hc <- private$cols$popuid
      ac <- private$cols$rel_abund

      # discrete rank colouring: one colour per haplotype rank (interpolating a ramp
      # from the palette when there are more ranks than colours), plus (if marked) a
      # separate category for invariant targets.
      if (!is.null(rank)) {
        has_inv <- "is_invariant" %in% names(prep_data) &&
          any(prep_data[["is_invariant"]], na.rm = TRUE)
        n_rank <- max(prep_data[["popid"]], na.rm = TRUE)
        if (n_rank > length(rank$palette)) {
          warning("More haplotype ranks (", n_rank, ") than palette colours (",
                  length(rank$palette),
                  "); interpolating a colour ramp from the palette.", call. = FALSE)
        }
        rank$palette <- .expand_rank_palette(rank$palette, n_rank)
        rank_chr <- as.character(prep_data[["popid"]])
        if (has_inv) rank_chr[which(prep_data[["is_invariant"]])] <- "invariant"
        lev <- c(as.character(seq_len(n_rank)), if (has_inv) "invariant")
        prep_data[[".rank"]] <- factor(rank_chr, levels = lev)
        fill_col <- ".rank"
      }

      # extra (non-standard) aesthetics carried through so ggplotly() can surface
      # them in hover text. Dynamic aesthetic names require !!! splicing in aes().
      hover <- rlang::set_names(
        list(rlang::sym(sc), rlang::sym(hc), rlang::sym(tc), rlang::sym(ac)),
        c(sc, hc, tc, ac)
      )
      mapping <- ggplot2::aes(
        xmin = as.numeric(.data[[tc]]) - 0.5,
        xmax = as.numeric(.data[[tc]]) + 0.5,
        ymin = as.numeric(.data[[sc]]) + .data[["fracModCumSum"]] - 0.5,
        ymax = as.numeric(.data[[sc]]) + .data[["fracModCumSum"]] +
          .data[["relAbundCol_mod"]] - 0.5,
        fill = .data[[fill_col]],
        !!!hover
      )

      p <- ggplot2::ggplot(prep_data) +
        ggplot2::geom_rect(mapping = mapping, color = "black") +
        .rainbow_theme()

      if (!is.null(rank)) {
        vals <- c(
          stats::setNames(rank$palette, as.character(seq_along(rank$palette))),
          invariant = rank$invariant_color
        )
        p <- p + ggplot2::scale_fill_manual(name = rank$title, values = vals,
                                            na.value = "grey80")
      } else {
        p <- p + ggplot2::guides(fill = "none") +
          if (style == "rainbow") {
            ggplot2::scale_fill_gradientn(colours = colors)
          } else {
            ggplot2::scale_fill_identity()
          }
      }

      target_levels <- levels(dplyr::pull(prep_data, dplyr::all_of(tc)))
      sample_levels <- levels(dplyr::pull(prep_data, dplyr::all_of(sc)))

      # x-axis: always expand = 0 (tight margins); drop ticks/labels when off.
      if (x_axis_labels) {
        p <- p +
          ggplot2::scale_x_continuous(
            labels = target_levels,
            breaks = seq_along(target_levels),
            expand = c(0, 0)
          ) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(
            family = "mono", angle = -90, hjust = 0
          ))
      } else {
        p <- p +
          ggplot2::scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank())
      }

      # y-axis: same treatment.
      if (y_axis_labels) {
        p <- p +
          ggplot2::scale_y_continuous(
            labels = sample_levels,
            breaks = seq_along(sample_levels),
            expand = c(0, 0)
          ) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(family = "mono"))
      } else {
        p <- p +
          ggplot2::scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
          ggplot2::theme(axis.title.y = ggplot2::element_blank())
      }
      p
    }
  )
)

#' Create a HaplotypeRainbow
#'
#' Convenience wrapper around `HaplotypeRainbow$new()`.
#'
#' @param data A data frame / tibble with one row per (sample, target, haplotype).
#' @param sample_col Name of the sample identifier column.
#' @param target_col Name of the target / locus column.
#' @param popuid_col Name of the within-target haplotype identifier column.
#' @param rel_abund_col Name of the relative counts column.
#' @return A [HaplotypeRainbow] object.
#' @examples
#' \dontrun{
#' rb <- haplotype_rainbow(my_allele_table)
#' }
#' @export
haplotype_rainbow <- function(data,
                              sample_col    = "library_sample_name",
                              target_col    = "target_name",
                              popuid_col    = "seq",
                              rel_abund_col = "reads") {
  HaplotypeRainbow$new(data, sample_col, target_col, popuid_col, rel_abund_col)
}
