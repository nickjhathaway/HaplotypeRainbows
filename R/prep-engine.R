# Internal prep engine ---------------------------------------------------------
#
# The three historical prep functions (prepForRainbow, prepForRainbowArrangedByFrac
# and prepForRainbowShade) were ~90% identical. They are unified here into a single
# engine. The user's four mapped columns are pulled out of the input (every other
# column is dropped) and given canonical internal names -- sample, target, hapid and
# rel_abund -- so the engine, and every method that later reads the prepped table,
# works against a small, fixed, collision-free schema regardless of the user's column
# names. `rel_abund` keeps the RAW counts; the within-sample/target fraction is added
# as `within_sample_freq` (with `total_abund` the per-(sample, target) denominator).
# All derived columns are snake_case.

# Select ONLY the four mapped columns and give them canonical internal names.
# transmute() drops every other input column, which keeps the internal frame small and
# -- crucially -- avoids name collisions when a user column happens to share a canonical
# name (e.g. a stray SeekDeep `p_name` left in the data when the target was remapped to
# a different column).
.canonicalize <- function(data, cols) {
  dplyr::transmute(
    data,
    sample    = .data[[cols$sample]],
    target    = .data[[cols$target]],
    hapid     = .data[[cols$popuid]],
    rel_abund = .data[[cols$rel_abund]]
  )
}

# Shared core: collapse to (sample, target, haplotype), keep the raw counts in
# `rel_abund`, add the per-(sample, target) `total_abund` and the `within_sample_freq`
# fraction, and compute the cumulative-sum bar geometry (from the fraction). When
# `arrange_by_frac` is TRUE the rows are ordered by ascending fraction before the
# cumulative sums (this is the only difference between the "population" and "frac"
# preps).
.prep_core <- function(data, bar_height, arrange_by_frac) {
  out <- data %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(n_targets = length(unique(target))) %>%
    dplyr::group_by() %>%
    dplyr::mutate(sample = as.character(sample)) %>%
    dplyr::mutate(sample = factor(sample)) %>%
    dplyr::group_by(sample) %>%
    dplyr::arrange(hapid) %>%
    dplyr::group_by(sample, target, hapid) %>%
    dplyr::summarise(rel_abund = sum(rel_abund), .groups = "drop_last") %>%
    dplyr::group_by(sample, target) %>%
    dplyr::mutate(total_abund = sum(rel_abund)) %>%
    dplyr::mutate(within_sample_freq = rel_abund / total_abund)

  if (arrange_by_frac) {
    out <- out %>% dplyr::arrange(within_sample_freq)
  }

  out <- out %>%
    dplyr::group_by(sample, target, hapid) %>%
    dplyr::mutate(sample_coi = length(unique(hapid))) %>%
    dplyr::group_by(sample, target) %>%
    dplyr::mutate(
      within_sample_freq_mod = within_sample_freq * bar_height,
      freq_cumsum          = cumsum(within_sample_freq) - within_sample_freq,
      freq_mod_cumsum      = cumsum(within_sample_freq_mod) - within_sample_freq_mod,
      fake_freq            = 1 / unique(sample_coi),
      fake_freq_mod        = fake_freq * bar_height,
      fake_freq_cumsum     = cumsum(fake_freq) - fake_freq,
      fake_freq_mod_cumsum = cumsum(fake_freq_mod) - fake_freq_mod
    )
  out
}

# Rank haplotypes within each target by how many samples carry them.
.prep_popname <- function(core) {
  core %>%
    dplyr::select(sample, target, hapid) %>%
    unique() %>%
    dplyr::group_by(target, hapid) %>%
    dplyr::summarise(samp_n = dplyr::n(), .groups = "drop_last") %>%
    dplyr::arrange(target, dplyr::desc(samp_n)) %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(pop_id = dplyr::row_number()) %>%
    dplyr::mutate(max_pop_id = max(pop_id))
}

# Population / frac colouring tail: rotate a hue offset across targets so that colours
# cycle with `color_period`, producing the rainbow.
.prep_rainbow_colors <- function(core, min_pop_size, color_period) {
  filt <- core %>%
    dplyr::left_join(.prep_popname(core),
                     by = c("target", "hapid")) %>%
    dplyr::filter(max_pop_id >= min_pop_size) %>%
    dplyr::group_by() %>%
    dplyr::mutate(target = factor(target))

  target_levels <- levels(dplyr::pull(filt, target))
  target_to_hue <- tibble::tibble(
    target = factor(target_levels, levels = target_levels),
    hue_mod = (seq_along(target_levels) - 1) %% color_period + 1
  )

  filt %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(pop_id_frac = (pop_id - 1) / max_pop_id) %>%
    dplyr::left_join(target_to_hue, by = "target") %>%
    dplyr::mutate(
      pop_id_perc           = 100 * pop_id_frac,
      pop_id_frac_reg_color = round(abs((pop_id_perc + (hue_mod / color_period) * 100) %% 200 - 0.0001) %% 100),
      pop_id_perc_log       = log((pop_id_frac * 99) + 1, base = 100) * 100,
      pop_id_frac_log_color = round(abs((pop_id_perc_log + (hue_mod / color_period) * 100) %% 200 - 0.0001) %% 100)
    ) %>%
    dplyr::group_by()
}

# Shade colouring tail: assign each target a base colour and shade haplotypes within it.
.prep_shade_colors <- function(core, min_pop_size, base_colors) {
  filt <- core %>%
    dplyr::left_join(.prep_popname(core),
                     by = c("target", "hapid")) %>%
    dplyr::filter(max_pop_id >= min_pop_size) %>%
    dplyr::group_by() %>%
    dplyr::mutate(target = factor(target))

  pop_colors <- core %>%
    dplyr::left_join(.prep_popname(core),
                     by = c("target", "hapid")) %>%
    dplyr::group_by(target, hapid) %>%
    dplyr::count() %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(total = sum(n)) %>%
    dplyr::arrange(target, n) %>%
    dplyr::mutate(freq = n / total) %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(
      p_uniq_haps   = dplyr::n(),
      h_id          = dplyr::row_number(),
      h_id_freq     = h_id / p_uniq_haps,
      h_id_freq_mod = h_id_freq * 0.75 + 0.25,
      cum_freq      = cumsum(freq),
      mod_cum_freq  = cum_freq * 0.75 + 0.25,
      target        = factor(target),
      p_color_id    = (as.numeric(target) %% length(base_colors)) + 1,
      p_hue         = p_color_id / length(base_colors),
      p_base_color  = base_colors[p_color_id],
      h_color             = grDevices::hsv(p_hue, alpha = cum_freq),
      h_color_mod         = scales::alpha(p_base_color, alpha = mod_cum_freq),
      h_color_by_freq     = scales::alpha(p_base_color, alpha = h_id_freq),
      h_color_by_freq_mod = scales::alpha(p_base_color, alpha = h_id_freq_mod)
    )

  filt %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(pop_id_frac = (pop_id - 1) / max_pop_id) %>%
    dplyr::left_join(pop_colors, by = c("target", "hapid")) %>%
    dplyr::group_by()
}

# Top-level engine used by the HaplotypeRainbow class. Returns the prepped table with
# canonical column names (the user's extra columns are intentionally not carried).
.prep_engine <- function(data, cols,
                         sort = c("population_rank", "within_sample_freq", "shade"),
                         min_pop_size = 1, color_period = 11, bar_height = 0.80,
                         base_colors = c("#e41a1c", "#377eb8", "#4daf4a",
                                         "#984ea3", "#ff7f00", "#ffff33")) {
  sort <- match.arg(sort)
  data <- .canonicalize(data, cols)
  core <- .prep_core(data, bar_height,
                     arrange_by_frac = (sort == "within_sample_freq"))

  switch(
    sort,
    population_rank    = .prep_rainbow_colors(core, min_pop_size, color_period),
    within_sample_freq = .prep_rainbow_colors(core, min_pop_size, color_period),
    shade              = .prep_shade_colors(core, min_pop_size, base_colors)
  )
}
