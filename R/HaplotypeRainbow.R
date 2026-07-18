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

    #' @description Prep the data for plotting. Stores the prepped table internally.
    #' @param sort One of "population" (haplotypes ordered by population rank),
    #'   "frac" (ordered by within-sample fraction) or "shade" (colour by shading a
    #'   per-target base colour instead of a rainbow).
    #' @param min_pop_size Drop targets with fewer than this many unique haplotypes.
    #' @param color_period Number of hue steps to rotate across targets (rainbow modes).
    #' @param bar_height Height of the full stacked bar per sample; < 1 leaves a gap.
    #' @param base_colors Base colours to shade from (shade mode only).
    #' @return The object, invisibly (chainable).
    prep = function(sort = c("population", "frac", "shade"),
                    min_pop_size = 1, color_period = 11, bar_height = 0.80,
                    base_colors = c("#e41a1c", "#377eb8", "#4daf4a",
                                    "#984ea3", "#ff7f00", "#ffff33")) {
      sort <- match.arg(sort)
      private$prepped <- .prep_engine(
        private$data, private$cols, sort = sort,
        min_pop_size = min_pop_size, color_period = color_period,
        bar_height = bar_height, base_colors = base_colors
      )
      private$sort_mode <- sort
      private$color_period <- color_period
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

    #' @description Build the rainbow plot.
    #' @param style "rainbow" (gradient over rotating hues) or "shade" (identity colours
    #'   from shade prep). Defaults to "shade" when the data was prepped with
    #'   `sort = "shade"`, otherwise "rainbow".
    #' @param colors Colour ramp for rainbow style. Length should match `color_period`.
    #' @param color_col Column to map to fill for rainbow style
    #'   ("popidFracLogColor" or "popidFracRegColor").
    #' @param shade_col Identity-colour column for shade style.
    #' @param axis_labels Add target names on the x-axis and sample names on the y-axis.
    #' @return A [ggplot2::ggplot] object.
    plot = function(style = NULL,
                    colors = RColorBrewer::brewer.pal(11, "Spectral"),
                    color_col = "popidFracLogColor",
                    shade_col = "h_color_byFreq_mod",
                    axis_labels = TRUE) {
      private$require_prepped()
      if (is.null(style)) {
        style <- if (identical(private$sort_mode, "shade")) "shade" else "rainbow"
      }
      style <- match.arg(style, c("rainbow", "shade"))
      if (style == "rainbow" && identical(private$sort_mode, "shade")) {
        stop("style = 'rainbow' needs data prepped with sort = 'population' or 'frac'.",
             call. = FALSE)
      }
      if (style == "shade" && !identical(private$sort_mode, "shade")) {
        stop("style = 'shade' needs data prepped with sort = 'shade'.", call. = FALSE)
      }
      fill_col <- if (style == "rainbow") color_col else shade_col
      private$build_plot(style, fill_col, colors, axis_labels)
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

    require_prepped = function() {
      if (is.null(private$prepped)) {
        stop("Call $prep() before this operation.", call. = FALSE)
      }
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

    build_plot = function(style, fill_col, colors, axis_labels) {
      prep_data <- private$prepped
      sc <- private$cols$sample
      tc <- private$cols$target
      hc <- private$cols$popuid
      ac <- private$cols$rel_abund

      y_axis <- .sample_axis(prep_data, sc)

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

      if (axis_labels) {
        target_levels <- levels(dplyr::pull(prep_data, dplyr::all_of(tc)))
        sample_levels <- levels(dplyr::pull(prep_data, dplyr::all_of(sc)))
        p <- p +
          ggplot2::scale_x_continuous(
            labels = target_levels,
            breaks = seq_along(target_levels),
            expand = c(0, 0)
          ) +
          ggplot2::scale_y_continuous(
            labels = sample_levels,
            breaks = seq_along(sample_levels),
            expand = c(0, 0)
          ) +
          ggplot2::theme(
            axis.text.y = ggplot2::element_text(family = "mono"),
            axis.text.x = ggplot2::element_text(family = "mono", angle = -90, hjust = 0)
          )
      } else {
        p <- p +
          ggplot2::scale_y_continuous(breaks = y_axis$breaks, labels = y_axis$labels) +
          ggplot2::theme(axis.text.x = ggplot2::element_blank())
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
