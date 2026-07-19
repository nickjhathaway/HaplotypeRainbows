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
    #' @return The object, invisibly (chainable).
    prep = function(sort = c("population_rank", "within_sample_freq"),
                    min_pop_size = 1, color_period = 11, bar_height = 0.80) {
      sort <- match.arg(sort)
      private$prepped <- .prep_engine(
        private$data, private$cols, sort = sort,
        min_pop_size = min_pop_size, color_period = color_period,
        bar_height = bar_height
      )
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
    #' @param x_axis_labels Show target names on the x-axis.
    #' @param y_axis_labels Show sample names on the y-axis.
    #' @return A [ggplot2::ggplot] object.
    plot = function(style = NULL,
                    colors = RColorBrewer::brewer.pal(11, "Spectral"),
                    color_col = "popidFracLogColor",
                    shade_col = "h_color_byFreq_mod",
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
      fill_col <- if (style == "rainbow") color_col else shade_col
      private$build_plot(style, fill_col, colors, x_axis_labels, y_axis_labels)
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

    build_plot = function(style, fill_col, colors, x_axis_labels, y_axis_labels) {
      prep_data <- private$prepped
      sc <- private$cols$sample
      tc <- private$cols$target
      hc <- private$cols$popuid
      ac <- private$cols$rel_abund

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
        .rainbow_theme() +
        ggplot2::guides(fill = "none")

      p <- p + if (style == "rainbow") {
        ggplot2::scale_fill_gradientn(colours = colors)
      } else {
        ggplot2::scale_fill_identity()
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
