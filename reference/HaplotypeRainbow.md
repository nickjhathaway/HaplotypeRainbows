# HaplotypeRainbow

An R6 class for building haplotype "rainbow" plots. The class carries
the mapping of your four data columns (sample, target, haplotype
identifier and relative counts) so it only has to be supplied once.
Methods that transform state (`prep()`, `sort_by_clustering()`,
`set_sample_order()`, `sort_alphabetical()`) mutate the object and
return it invisibly, so they can be chained;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

Column defaults follow the Portable Microhaplotype Object (PMO)
convention: `library_sample_name`, `target_name`, `seq` and `reads`.

## Methods

### Public methods

- [`HaplotypeRainbow$new()`](#method-HaplotypeRainbow-initialize)

- [`HaplotypeRainbow$prep()`](#method-HaplotypeRainbow-prep)

- [`HaplotypeRainbow$prep_shade()`](#method-HaplotypeRainbow-prep_shade)

- [`HaplotypeRainbow$sort_by_clustering()`](#method-HaplotypeRainbow-sort_by_clustering)

- [`HaplotypeRainbow$get_hclust()`](#method-HaplotypeRainbow-get_hclust)

- [`HaplotypeRainbow$get_dendrogram()`](#method-HaplotypeRainbow-get_dendrogram)

- [`HaplotypeRainbow$cluster_groups()`](#method-HaplotypeRainbow-cluster_groups)

- [`HaplotypeRainbow$set_sample_order()`](#method-HaplotypeRainbow-set_sample_order)

- [`HaplotypeRainbow$sort_alphabetical()`](#method-HaplotypeRainbow-sort_alphabetical)

- [`HaplotypeRainbow$set_sample_meta()`](#method-HaplotypeRainbow-set_sample_meta)

- [`HaplotypeRainbow$set_target_meta()`](#method-HaplotypeRainbow-set_target_meta)

- [`HaplotypeRainbow$get_sample_meta()`](#method-HaplotypeRainbow-get_sample_meta)

- [`HaplotypeRainbow$get_target_meta()`](#method-HaplotypeRainbow-get_target_meta)

- [`HaplotypeRainbow$sort_samples_by_meta()`](#method-HaplotypeRainbow-sort_samples_by_meta)

- [`HaplotypeRainbow$sort_targets_by_meta()`](#method-HaplotypeRainbow-sort_targets_by_meta)

- [`HaplotypeRainbow$set_target_order()`](#method-HaplotypeRainbow-set_target_order)

- [`HaplotypeRainbow$plot()`](#method-HaplotypeRainbow-plot)

- [`HaplotypeRainbow$get_rank_colors()`](#method-HaplotypeRainbow-get_rank_colors)

- [`HaplotypeRainbow$add_sample_metadata()`](#method-HaplotypeRainbow-add_sample_metadata)

- [`HaplotypeRainbow$add_target_annotation()`](#method-HaplotypeRainbow-add_target_annotation)

- [`HaplotypeRainbow$dims()`](#method-HaplotypeRainbow-dims)

- [`HaplotypeRainbow$save_pdf()`](#method-HaplotypeRainbow-save_pdf)

- [`HaplotypeRainbow$drop_legends()`](#method-HaplotypeRainbow-drop_legends)

- [`HaplotypeRainbow$extract_legend()`](#method-HaplotypeRainbow-extract_legend)

- [`HaplotypeRainbow$save_legend_pdf()`](#method-HaplotypeRainbow-save_legend_pdf)

- [`HaplotypeRainbow$export_groups_pdf()`](#method-HaplotypeRainbow-export_groups_pdf)

- [`HaplotypeRainbow$add_cluster_gaps()`](#method-HaplotypeRainbow-add_cluster_gaps)

- [`HaplotypeRainbow$add_cluster_bands()`](#method-HaplotypeRainbow-add_cluster_bands)

- [`HaplotypeRainbow$get_prepped()`](#method-HaplotypeRainbow-get_prepped)

- [`HaplotypeRainbow$column_map()`](#method-HaplotypeRainbow-column_map)

- [`HaplotypeRainbow$print()`](#method-HaplotypeRainbow-print)

- [`HaplotypeRainbow$clone()`](#method-HaplotypeRainbow-clone)

------------------------------------------------------------------------

### `HaplotypeRainbow$new()`

Create a new HaplotypeRainbow.

#### Usage

    HaplotypeRainbow$new(
      data,
      sample_col = "library_sample_name",
      target_col = "target_name",
      popuid_col = "seq",
      rel_abund_col = "reads"
    )

#### Arguments

- `data`:

  A data frame / tibble with one row per (sample, target, haplotype).

- `sample_col`:

  Name of the sample identifier column.

- `target_col`:

  Name of the target / locus column.

- `popuid_col`:

  Name of the within-target haplotype identifier column.

- `rel_abund_col`:

  Name of the relative counts column (raw counts are fine; they are
  normalised to within-sample/target fractions internally).

------------------------------------------------------------------------

### `HaplotypeRainbow$prep()`

Prep the data for a rainbow plot (haplotype colours rotate across
targets). Stores the prepped table internally.

#### Usage

    HaplotypeRainbow$prep(
      sort = c("population_rank", "within_sample_freq"),
      min_pop_size = 1,
      color_period = 11,
      bar_height = 0.8,
      mark_invariant = FALSE
    )

#### Arguments

- `sort`:

  Haplotype ordering within each cell: "population_rank" (order by
  population rank) or "within_sample_freq" (order by within-sample
  fraction).

- `min_pop_size`:

  Drop targets with fewer than this many unique haplotypes.

- `color_period`:

  Number of hue steps to rotate across targets.

- `bar_height`:

  Height of the full stacked bar per sample; \< 1 leaves a gap.

- `mark_invariant`:

  Flag single-haplotype (invariant) targets so that
  `plot(rank_colors = TRUE)` can colour them separately (black by
  default).

#### Returns

The object, invisibly (chainable).

------------------------------------------------------------------------

### `HaplotypeRainbow$prep_shade()`

Prep the data for a shade plot (haplotypes shaded from a per-target base
colour instead of a rotating rainbow). Stores the prepped table
internally.

#### Usage

    HaplotypeRainbow$prep_shade(
      min_pop_size = 1,
      bar_height = 0.8,
      base_colors = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")
    )

#### Arguments

- `min_pop_size`:

  Drop targets with fewer than this many unique haplotypes.

- `bar_height`:

  Height of the full stacked bar per sample; \< 1 leaves a gap.

- `base_colors`:

  Base colours to generate the per-target shades from.

#### Returns

The object, invisibly (chainable).

------------------------------------------------------------------------

### `HaplotypeRainbow$sort_by_clustering()`

Reorder samples so that similar samples (by haplotype sharing) sit next
to each other, via hierarchical clustering. The fitted `hclust` object
is stored (see `get_hclust()` / `get_dendrogram()` /
`cluster_groups()`).

#### Usage

    HaplotypeRainbow$sort_by_clustering(
      by_major_allele = FALSE,
      coverage_cutoff = 0.8,
      targets = NULL,
      dist_method = "euclidean",
      hclust_method = "ward.D2"
    )

#### Arguments

- `by_major_allele`:

  Cluster on only the major haplotype per (sample, target).

- `coverage_cutoff`:

  Keep only targets present in at least this fraction of samples for the
  clustering (ignored when `targets` is supplied).

- `targets`:

  Optional character vector of target names to cluster on; when given,
  only these targets are used (instead of the coverage-based selection).

- `dist_method`:

  Distance method passed to
  [`stats::dist()`](https://rdrr.io/r/stats/dist.html).

- `hclust_method`:

  Linkage method passed to
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html).

#### Returns

The object, invisibly (chainable).

------------------------------------------------------------------------

### `HaplotypeRainbow$get_hclust()`

Return the `hclust` object from the last `sort_by_clustering()`.

#### Usage

    HaplotypeRainbow$get_hclust()

#### Returns

An [stats::hclust](https://rdrr.io/r/stats/hclust.html) object.

------------------------------------------------------------------------

### `HaplotypeRainbow$get_dendrogram()`

Return the clustering as a dendrogram.

#### Usage

    HaplotypeRainbow$get_dendrogram()

#### Returns

A [stats::dendrogram](https://rdrr.io/r/stats/dendrogram.html).

------------------------------------------------------------------------

### `HaplotypeRainbow$cluster_groups()`

Cut the stored clustering into groups (by number of groups `k` or height
`h`). Groups are relabelled 1, 2, ... in dendrogram order.

#### Usage

    HaplotypeRainbow$cluster_groups(k = NULL, h = NULL)

#### Arguments

- `k`:

  Desired number of groups.

- `h`:

  Height at which to cut the tree.

#### Returns

A tibble with the sample column and a `cluster` factor.

------------------------------------------------------------------------

### `HaplotypeRainbow$set_sample_order()`

Set the sample order explicitly.

#### Usage

    HaplotypeRainbow$set_sample_order(levels)

#### Arguments

- `levels`:

  Character vector of sample names in the desired order.

#### Returns

The object, invisibly (chainable).

------------------------------------------------------------------------

### `HaplotypeRainbow$sort_alphabetical()`

Order samples alphabetically (the default).

#### Usage

    HaplotypeRainbow$sort_alphabetical()

#### Returns

The object, invisibly (chainable).

------------------------------------------------------------------------

### `HaplotypeRainbow$set_sample_meta()`

Attach per-sample metadata to the object (used by the metadata sidebar
and by `sort_samples_by_meta()`).

#### Usage

    HaplotypeRainbow$set_sample_meta(meta, match_col, cols = NULL, add = FALSE)

#### Arguments

- `meta`:

  A data frame of sample metadata.

- `match_col`:

  Name of the column in `meta` that holds the sample identifiers.

- `cols`:

  Metadata columns to keep (default: all columns except `match_col`).

- `add`:

  If `FALSE` (default) replace any existing sample metadata; if `TRUE`
  merge these columns onto the existing sample metadata (joined on the
  samples).

#### Returns

The object, invisibly (chainable).

------------------------------------------------------------------------

### `HaplotypeRainbow$set_target_meta()`

Attach per-target metadata to the object (used by the target annotation
strip and by `sort_targets_by_meta()`).

#### Usage

    HaplotypeRainbow$set_target_meta(meta, match_col, cols = NULL, add = FALSE)

#### Arguments

- `meta`:

  A data frame of target metadata.

- `match_col`:

  Name of the column in `meta` that holds the target names.

- `cols`:

  Metadata columns to keep (default: all columns except `match_col`).

- `add`:

  If `FALSE` (default) replace any existing target metadata; if `TRUE`
  merge these columns onto the existing target metadata (joined on the
  targets).

#### Returns

The object, invisibly (chainable).

------------------------------------------------------------------------

### `HaplotypeRainbow$get_sample_meta()`

Return the stored sample metadata (or `NULL`).

#### Usage

    HaplotypeRainbow$get_sample_meta()

#### Returns

A data frame.

------------------------------------------------------------------------

### `HaplotypeRainbow$get_target_meta()`

Return the stored target metadata (or `NULL`).

#### Usage

    HaplotypeRainbow$get_target_meta()

#### Returns

A data frame.

------------------------------------------------------------------------

### `HaplotypeRainbow$sort_samples_by_meta()`

Order samples by one or more sample-metadata columns. Samples with no
metadata are placed at the end.

#### Usage

    HaplotypeRainbow$sort_samples_by_meta(cols, desc = FALSE)

#### Arguments

- `cols`:

  Character vector of sample-metadata column names to order by.

- `desc`:

  Sort descending instead of ascending.

#### Returns

The object, invisibly (chainable).

------------------------------------------------------------------------

### `HaplotypeRainbow$sort_targets_by_meta()`

Order targets by one or more target-metadata columns. Targets with no
metadata are placed at the end.

#### Usage

    HaplotypeRainbow$sort_targets_by_meta(cols, desc = FALSE)

#### Arguments

- `cols`:

  Character vector of target-metadata column names to order by.

- `desc`:

  Sort descending instead of ascending.

#### Returns

The object, invisibly (chainable).

------------------------------------------------------------------------

### `HaplotypeRainbow$set_target_order()`

Set the target (column) order explicitly.

#### Usage

    HaplotypeRainbow$set_target_order(levels)

#### Arguments

- `levels`:

  Character vector of target names in the desired order.

#### Returns

The object, invisibly (chainable).

------------------------------------------------------------------------

### `HaplotypeRainbow$plot()`

Build the rainbow plot.

#### Usage

    HaplotypeRainbow$plot(
      style = NULL,
      colors = RColorBrewer::brewer.pal(11, "Spectral"),
      color_col = "pop_id_frac_log_color",
      shade_col = "h_color_by_freq_mod",
      rank_colors = FALSE,
      rank_palette = colorPalette_12,
      invariant_color = "#000000",
      rank_legend_title = "Microhaplotype Rank",
      x_axis_labels = TRUE,
      y_axis_labels = TRUE
    )

#### Arguments

- `style`:

  "rainbow" (gradient over rotating hues) or "shade" (identity colours
  from shade prep). Defaults to "shade" when the data was prepped with
  `sort = "shade"`, otherwise "rainbow".

- `colors`:

  Colour ramp for rainbow style. Length should match `color_period`.

- `color_col`:

  Column to map to fill for rainbow style ("pop_id_frac_log_color" or
  "pop_id_frac_reg_color").

- `shade_col`:

  Identity-colour column for shade style.

- `rank_colors`:

  Colour haplotypes by a discrete rank (with a legend) instead of the
  continuous rainbow gradient. When there are more ranks than palette
  colours, a colour ramp is interpolated from the palette (with a
  warning); invariant targets (see `prep(mark_invariant = TRUE)`) get
  `invariant_color`. Replaces the old harvest-fills-from-`ggplot_build`
  recipe.

- `rank_palette`:

  Colours for haplotype ranks 1, 2, ... (rank colouring); if there are
  more ranks than colours, a ramp is interpolated from these.

- `invariant_color`:

  Colour for invariant (single-haplotype) targets.

- `rank_legend_title`:

  Legend title used with `rank_colors`.

- `x_axis_labels`:

  Show target names on the x-axis.

- `y_axis_labels`:

  Show sample names on the y-axis.

#### Returns

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

------------------------------------------------------------------------

### `HaplotypeRainbow$get_rank_colors()`

Return the discrete rank -\> colour map used by
`plot(rank_colors = TRUE)` (haplotype ranks plus the invariant colour).
When there are more ranks than palette colours a ramp is interpolated.

#### Usage

    HaplotypeRainbow$get_rank_colors(
      n_ranks = NULL,
      rank_palette = colorPalette_12,
      invariant_color = "#000000"
    )

#### Arguments

- `n_ranks`:

  Number of ranks to produce colours for (default: the maximum haplotype
  rank in the prepped data, or the palette length if not prepped).

- `rank_palette`:

  Colours for ranks 1, 2, ...

- `invariant_color`:

  Colour for invariant targets.

#### Returns

A named character vector.

------------------------------------------------------------------------

### `HaplotypeRainbow$add_sample_metadata()`

Add a per-sample metadata sidebar (coloured bands) to a rainbow plot,
using the metadata attached with `set_sample_meta()`.

#### Usage

    HaplotypeRainbow$add_sample_metadata(
      p,
      cols = NULL,
      side = "left",
      width = 1,
      height = 1,
      gap = 0,
      plot_gap = 1,
      colors = NULL,
      na_color = "grey80",
      legend = TRUE,
      legend_ncol = NULL,
      legend_nrow = NULL,
      level_order = NULL,
      labels = TRUE,
      target_labels = TRUE,
      border = "black"
    )

#### Arguments

- `p`:

  A ggplot returned by
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- `cols`:

  Metadata columns to draw (default: all in the sample-metadata slot),
  one band each, in order outward from the plot.

- `side`:

  "left" or "right" of the rainbow.

- `width`:

  Band width in cell units (1 = one rainbow cell wide).

- `height`:

  Band height as a fraction of a cell (1 = full cell, centred).

- `gap`:

  Gap between bands in cell units (0 = flush).

- `plot_gap`:

  Gap in cell units between the rainbow cells and the first band (0 =
  flush against the plot).

- `colors`:

  Optional named list (per column) of value -\> colour to override the
  auto-assigned colours.

- `na_color`:

  Fill for samples missing that metadata value.

- `legend`:

  Draw a colour legend for each band.

- `legend_ncol, legend_nrow`:

  Number of columns / rows for the band legends, to keep tall legends on
  the page. A scalar applies to every band; a named vector (e.g.
  `c(country = 3)`) sets it per column.

- `level_order`:

  Optional named list giving the legend/level order for a column (e.g.
  `list(country = c("Uganda", "Kenya"))`); unlisted values are appended
  in sorted order.

- `labels`:

  Draw the band (column-name) labels on the x-axis. Independent of
  `target_labels`, so band labels can stay on when target names are off.

- `target_labels`:

  Keep the target names on the x-axis.

- `border`:

  Border colour for the band rectangles.

#### Returns

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

------------------------------------------------------------------------

### `HaplotypeRainbow$add_target_annotation()`

Add a per-target annotation strip (coloured bands) above or below a
rainbow plot, using the metadata attached with `set_target_meta()`.

#### Usage

    HaplotypeRainbow$add_target_annotation(
      p,
      cols = NULL,
      position = "top",
      width = 1,
      height = 1,
      gap = 0,
      plot_gap = 1,
      colors = NULL,
      na_color = "grey80",
      legend = TRUE,
      legend_ncol = NULL,
      legend_nrow = NULL,
      level_order = NULL,
      labels = TRUE,
      sample_labels = TRUE,
      border = "black"
    )

#### Arguments

- `p`:

  A ggplot returned by
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- `cols`:

  Metadata columns to draw (default: all in the target-metadata slot),
  one band each, stacked outward from the plot.

- `position`:

  "top" or "bottom" of the rainbow.

- `width`:

  Band width as a fraction of a cell (1 = full cell, centred).

- `height`:

  Band thickness in cell units (1 = one rainbow cell tall).

- `gap`:

  Gap between bands in cell units (0 = flush).

- `plot_gap`:

  Gap in cell units between the rainbow cells and the first band (0 =
  flush against the plot).

- `colors`:

  Optional named list (per column) of value -\> colour to override the
  auto-assigned colours.

- `na_color`:

  Fill for targets missing that metadata value.

- `legend`:

  Draw a colour legend for each band.

- `legend_ncol, legend_nrow`:

  Number of columns / rows for the band legends. A scalar applies to
  every band; a named vector sets it per column.

- `level_order`:

  Optional named list giving the legend/level order for a column;
  unlisted values are appended in sorted order.

- `labels`:

  Draw the band (column-name) labels on the y-axis. Independent of
  `sample_labels`.

- `sample_labels`:

  Keep the sample names on the y-axis.

- `border`:

  Border colour for the band rectangles.

#### Returns

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

------------------------------------------------------------------------

### `HaplotypeRainbow$dims()`

Suggested figure dimensions (inches) for the current data, scaling width
with the number of targets and height with the number of samples.

#### Usage

    HaplotypeRainbow$dims(
      p = NULL,
      cell_width = 0.3,
      cell_height = 0.3,
      min_width = 6,
      min_height = 6,
      extra_width = 0,
      extra_height = 0,
      size_include_legend = FALSE
    )

#### Arguments

- `p`:

  A ggplot; only needed when `size_include_legend = TRUE`.

- `cell_width`:

  Inches per target column.

- `cell_height`:

  Inches per sample row.

- `min_width`:

  Minimum width.

- `min_height`:

  Minimum height.

- `extra_width`:

  Padding added to width (e.g. for a sidebar).

- `extra_height`:

  Padding added to height (e.g. for a target strip).

- `size_include_legend`:

  If `TRUE`, enlarge the size to fit the plot's (bottom) legend as well:
  widen to at least the legend width and add the legend height. Requires
  `p`.

#### Returns

A list with `width` and `height` (inches).

------------------------------------------------------------------------

### `HaplotypeRainbow$save_pdf()`

Save a plot to PDF, sized automatically from the data when `width` /
`height` are not supplied. Uses `cairo_pdf` by default (which keeps
hyphens as hyphens rather than converting them to en/em dashes); pass
`device = "pdf"` for the base device (e.g. on Windows, where the cairo
device can misbehave).

#### Usage

    HaplotypeRainbow$save_pdf(
      p,
      file,
      width = NULL,
      height = NULL,
      device = c("cairo", "pdf"),
      size_include_legend = FALSE,
      bg = "transparent",
      ...
    )

#### Arguments

- `p`:

  A ggplot (e.g. from
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) or
  `add_sample_metadata()`).

- `file`:

  Output file path.

- `width`:

  Width in inches (default from `dims()`).

- `height`:

  Height in inches (default from `dims()`).

- `device`:

  "cairo" (default) or "pdf"; cairo falls back to pdf() when the
  platform lacks cairo support.

- `size_include_legend`:

  If `TRUE` (and sizing automatically), enlarge the figure so the plot's
  legend fits too. Default `FALSE`.

- `bg`:

  Background colour of the device; defaults to "transparent".

- `...`:

  Passed to the graphics device.

#### Returns

The file path, invisibly.

------------------------------------------------------------------------

### `HaplotypeRainbow$drop_legends()`

Return a copy of a plot with all legends removed. Combine with
`save_legend_pdf()` to export the plot and its legend separately (e.g.
to place the legend elsewhere in a figure).

#### Usage

    HaplotypeRainbow$drop_legends(p)

#### Arguments

- `p`:

  A ggplot.

#### Returns

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object with `legend.position = "none"`.

------------------------------------------------------------------------

### `HaplotypeRainbow$extract_legend()`

Extract the legend (guide-box) of a plot as a grob, for composing it
elsewhere (e.g. with patchwork / cowplot / grid).

#### Usage

    HaplotypeRainbow$extract_legend(p)

#### Arguments

- `p`:

  A ggplot.

#### Returns

A grid grob (gtable).

------------------------------------------------------------------------

### `HaplotypeRainbow$save_legend_pdf()`

Save just the legend of a plot to PDF (so it can be combined in post
with the plot exported via `save_pdf(drop_legends(p), ...)`). When
`width` / `height` are not supplied, the legend's natural size is
estimated so it is not clipped.

#### Usage

    HaplotypeRainbow$save_legend_pdf(
      p,
      file,
      width = NULL,
      height = NULL,
      margin = 0.2,
      device = c("cairo", "pdf"),
      bg = "transparent",
      ...
    )

#### Arguments

- `p`:

  A ggplot.

- `file`:

  Output file path.

- `width`:

  Width in inches (default: estimated from the legend).

- `height`:

  Height in inches (default: estimated from the legend).

- `margin`:

  Inches of padding added around the estimated size.

- `device`:

  "cairo" (default) or "pdf"; cairo falls back to pdf() when the
  platform lacks cairo support.

- `bg`:

  Background colour of the device; defaults to "transparent".

- `...`:

  Passed to the graphics device.

#### Returns

The file path, invisibly.

------------------------------------------------------------------------

### `HaplotypeRainbow$export_groups_pdf()`

Export one plot per group (cluster or sample-metadata value) to PDF,
auto-sizing each plot to its own sample count. By default the per-group
pages are combined into a single multi-page PDF (with correctly-sized
pages) using the qpdf package.

#### Usage

    HaplotypeRainbow$export_groups_pdf(
      file,
      plot_fun,
      by = c("cluster", "meta"),
      k = NULL,
      h = NULL,
      meta_col = NULL,
      device = c("cairo", "pdf"),
      combine = TRUE,
      align_targets = TRUE,
      size_include_legend = FALSE,
      bg = "transparent",
      ...
    )

#### Arguments

- `file`:

  Output PDF path (or, with `combine = FALSE`, a template whose `.pdf`
  is suffixed with each group label).

- `plot_fun`:

  A function of one argument: it receives a per-group HaplotypeRainbow
  (with the prepped data filtered to that group) and returns a ggplot.
  e.g. `function(sub) sub$plot()`.

- `by`:

  Group by "cluster" (needs a prior `sort_by_clustering()` and `k`/`h`)
  or "meta" (needs `set_sample_meta()` and `meta_col`).

- `k, h`:

  Cluster cut (when `by = "cluster"`).

- `meta_col`:

  Sample-metadata column to group by (when `by = "meta"`).

- `device`:

  "cairo" (default) or "pdf"; cairo falls back to pdf() when the
  platform lacks cairo support.

- `combine`:

  Combine the per-group pages into one PDF (needs qpdf); if `FALSE`,
  write a separate file per group.

- `align_targets`:

  Pad sample-name labels to a common width (mono font) so the plot body
  / targets line up at the same position on every page. Default `TRUE`.

- `size_include_legend`:

  Passed to `dims()` for per-group sizing.

- `bg`:

  Device background colour (default "transparent").

- `...`:

  Passed to the graphics device.

#### Returns

The output path(s), invisibly.

------------------------------------------------------------------------

### `HaplotypeRainbow$add_cluster_gaps()`

Insert empty spacer rows between clusters, producing a real physical gap
between cluster blocks in the plot (the cells, sidebar and cluster bands
all shift together, since everything is positioned by the sample
factor). Requires a prior `sort_by_clustering()`. Call with `gap = 0` to
remove the spacers.

#### Usage

    HaplotypeRainbow$add_cluster_gaps(k = NULL, h = NULL, gap = 1)

#### Arguments

- `k`:

  Number of cluster groups (passed to `cluster_groups()`).

- `h`:

  Height at which to cut the tree (alternative to `k`).

- `gap`:

  Number of blank rows to insert between adjacent clusters.

#### Returns

The object, invisibly (chainable).

------------------------------------------------------------------------

### `HaplotypeRainbow$add_cluster_bands()`

Overlay alternating translucent background bands, one per cluster group,
spanning the plot width — a quick way to show clusters on a single plot.
Requires a prior `sort_by_clustering()`.

#### Usage

    HaplotypeRainbow$add_cluster_bands(
      p,
      k = NULL,
      h = NULL,
      colors = c("#00000025", "#AAAAAA25"),
      extend_left = 0,
      extend_right = 0,
      expand = 0,
      border = NA
    )

#### Arguments

- `p`:

  A ggplot returned by
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- `k`:

  Number of cluster groups (passed to `cluster_groups()`).

- `h`:

  Height at which to cut the tree (alternative to `k`).

- `colors`:

  Colours to alternate between (translucent by default so the rainbow
  shows through).

- `extend_left, extend_right`:

  Extend the bands (in cell units) to cover a sidebar or annotation
  strip.

- `expand`:

  Push the band edges this many cells past the plotting cells on both
  sides, widening the panel so the bands are visible as margin strips
  (the overlay alone can be hard to see over a dense cell grid).

- `border`:

  Border colour for the bands (default none). A visible colour draws a
  line around each cluster block.

#### Returns

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

------------------------------------------------------------------------

### `HaplotypeRainbow$get_prepped()`

Return the prepped data frame (or `NULL` if `prep()` not yet called).
The prepped table uses canonical column names (`sample`, `target`,
`hapid`, `rel_abund`); set `original_names = TRUE` to relabel those four
back to the column names supplied for the input data, so the result can
be rejoined to the source data (see `$column_map()`).

#### Usage

    HaplotypeRainbow$get_prepped(original_names = FALSE)

#### Arguments

- `original_names`:

  Relabel the canonical columns back to the user's originals.

#### Returns

A data frame.

------------------------------------------------------------------------

### `HaplotypeRainbow$column_map()`

Return the key mapping the canonical internal column names (`sample`,
`target`, `hapid`, `rel_abund`) used in the prepped table to the
original column names supplied for the input data. Use it to relabel
prepped data and reconnect it to the source input.

#### Usage

    HaplotypeRainbow$column_map()

#### Returns

A named character vector (names = canonical, values = original).

------------------------------------------------------------------------

### `HaplotypeRainbow$print()`

Print a short summary.

#### Usage

    HaplotypeRainbow$print(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### `HaplotypeRainbow$clone()`

The objects of this class are cloneable with this method.

#### Usage

    HaplotypeRainbow$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
rb <- HaplotypeRainbow$new(my_allele_table)
rb$prep(sort = "population")$sort_by_clustering()
rb$plot()
} # }
```
