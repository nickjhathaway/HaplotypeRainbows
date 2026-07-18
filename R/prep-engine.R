# Internal prep engine ---------------------------------------------------------
#
# The three historical prep functions (prepForRainbow, prepForRainbowArrangedByFrac
# and prepForRainbowShade) were ~90% identical. They are unified here into a single
# engine. To keep behaviour byte-identical to the original functions, the user's four
# columns are renamed to canonical internal names, the original dplyr pipeline is run
# verbatim against those canonical names, and the four columns are renamed back on the
# way out. Everything downstream (derived columns) is named the same regardless of the
# user's column mapping.

# Rename the four mapped columns to canonical internal names.
.canonicalize <- function(data, cols) {
  dplyr::rename(
    data,
    s_Sample       = dplyr::all_of(cols$sample),
    p_name         = dplyr::all_of(cols$target),
    h_popUID       = dplyr::all_of(cols$popuid),
    c_AveragedFrac = dplyr::all_of(cols$rel_abund)
  )
}

# Rename the four canonical columns back to the user's names.
.decanonicalize <- function(data, cols) {
  dplyr::rename(
    data,
    "{cols$sample}"    := "s_Sample",
    "{cols$target}"    := "p_name",
    "{cols$popuid}"    := "h_popUID",
    "{cols$rel_abund}" := "c_AveragedFrac"
  )
}

# Shared core: collapse to (sample, target, haplotype), normalise within-sample/target
# fraction, and compute the cumulative-sum bar geometry. When `arrange_by_frac` is TRUE
# the rows are ordered by ascending fraction before the cumulative sums (this is the only
# difference between the "population" and "frac" preps).
.prep_core <- function(data, bar_height, arrange_by_frac) {
  out <- data %>%
    dplyr::group_by(s_Sample) %>%
    dplyr::mutate(targetNumber = length(unique(p_name))) %>%
    dplyr::group_by() %>%
    dplyr::mutate(s_Sample = as.character(s_Sample)) %>%
    dplyr::mutate(s_Sample = factor(s_Sample)) %>%
    dplyr::group_by(s_Sample) %>%
    dplyr::arrange(h_popUID) %>%
    dplyr::group_by(s_Sample, p_name, h_popUID) %>%
    dplyr::summarise(c_AveragedFrac = sum(c_AveragedFrac), .groups = "drop_last") %>%
    dplyr::group_by(s_Sample, p_name) %>%
    dplyr::mutate(totalAbund = sum(c_AveragedFrac)) %>%
    dplyr::mutate(c_AveragedFrac = c_AveragedFrac / totalAbund)

  if (arrange_by_frac) {
    out <- out %>% dplyr::arrange(c_AveragedFrac)
  }

  out <- out %>%
    dplyr::group_by(s_Sample, p_name, h_popUID) %>%
    dplyr::mutate(s_COI = length(unique(h_popUID))) %>%
    dplyr::group_by(s_Sample, p_name) %>%
    dplyr::mutate(
      relAbundCol_mod   = c_AveragedFrac * bar_height,
      fracCumSum        = cumsum(c_AveragedFrac) - c_AveragedFrac,
      fracModCumSum     = cumsum(relAbundCol_mod) - relAbundCol_mod,
      fakeFrac          = 1 / unique(s_COI),
      fakeFracMod       = fakeFrac * bar_height,
      fakeFracCumSum    = cumsum(fakeFrac) - fakeFrac,
      fakeFracModCumSum = cumsum(fakeFracMod) - fakeFracMod
    )
  out
}

# Rank haplotypes within each target by how many samples carry them.
.prep_popname <- function(core) {
  core %>%
    dplyr::select(s_Sample, p_name, h_popUID) %>%
    unique() %>%
    dplyr::group_by(p_name, h_popUID) %>%
    dplyr::summarise(samp_n = dplyr::n(), .groups = "drop_last") %>%
    dplyr::arrange(p_name, dplyr::desc(samp_n)) %>%
    dplyr::group_by(p_name) %>%
    dplyr::mutate(popid = dplyr::row_number()) %>%
    dplyr::mutate(maxPopid = max(popid))
}

# Population / frac colouring tail: rotate a hue offset across targets so that colours
# cycle with `color_period`, producing the rainbow.
.prep_rainbow_colors <- function(core, min_pop_size, color_period) {
  filt <- core %>%
    dplyr::left_join(.prep_popname(core),
                     by = c("p_name", "h_popUID")) %>%
    dplyr::filter(maxPopid >= min_pop_size) %>%
    dplyr::group_by() %>%
    dplyr::mutate(p_name = factor(p_name))

  target_levels <- levels(dplyr::pull(filt, p_name))
  target_to_hue <- tibble::tibble(
    p_name = factor(target_levels, levels = target_levels),
    hueMod = (seq_along(target_levels) - 1) %% color_period + 1
  )

  filt %>%
    dplyr::group_by(p_name) %>%
    dplyr::mutate(popidFrac = (popid - 1) / maxPopid) %>%
    dplyr::left_join(target_to_hue, by = "p_name") %>%
    dplyr::mutate(
      popidPerc         = 100 * popidFrac,
      popidFracRegColor = round(abs((popidPerc + (hueMod / color_period) * 100) %% 200 - 0.0001) %% 100),
      popidPercLog      = log((popidFrac * 99) + 1, base = 100) * 100,
      popidFracLogColor = round(abs((popidPercLog + (hueMod / color_period) * 100) %% 200 - 0.0001) %% 100)
    ) %>%
    dplyr::group_by()
}

# Shade colouring tail: assign each target a base colour and shade haplotypes within it.
.prep_shade_colors <- function(core, min_pop_size, base_colors) {
  filt <- core %>%
    dplyr::left_join(.prep_popname(core),
                     by = c("p_name", "h_popUID")) %>%
    dplyr::filter(maxPopid >= min_pop_size) %>%
    dplyr::group_by() %>%
    dplyr::mutate(p_name = factor(p_name))

  pop_colors <- core %>%
    dplyr::left_join(.prep_popname(core),
                     by = c("p_name", "h_popUID")) %>%
    dplyr::group_by(p_name, h_popUID) %>%
    dplyr::count() %>%
    dplyr::group_by(p_name) %>%
    dplyr::mutate(total = sum(n)) %>%
    dplyr::arrange(p_name, n) %>%
    dplyr::mutate(freq = n / total) %>%
    dplyr::group_by(p_name) %>%
    dplyr::mutate(
      p_uniqHaps    = dplyr::n(),
      h_id          = dplyr::row_number(),
      h_id_freq     = h_id / p_uniqHaps,
      h_id_freq_mod = h_id_freq * 0.75 + 0.25,
      cumFreq       = cumsum(freq),
      modCumFreq    = cumFreq * 0.75 + 0.25,
      p_name        = factor(p_name),
      p_color_ID    = (as.numeric(p_name) %% length(base_colors)) + 1,
      p_hue         = p_color_ID / length(base_colors),
      p_baseColor   = base_colors[p_color_ID],
      h_color            = grDevices::hsv(p_hue, alpha = cumFreq),
      h_color_mod        = scales::alpha(p_baseColor, alpha = modCumFreq),
      h_color_byFreq     = scales::alpha(p_baseColor, alpha = h_id_freq),
      h_color_byFreq_mod = scales::alpha(p_baseColor, alpha = h_id_freq_mod)
    )

  filt %>%
    dplyr::group_by(p_name) %>%
    dplyr::mutate(popidFrac = (popid - 1) / maxPopid) %>%
    dplyr::left_join(pop_colors, by = c("p_name", "h_popUID")) %>%
    dplyr::group_by()
}

# Top-level engine used by the HaplotypeRainbow class.
.prep_engine <- function(data, cols, sort = c("population", "frac", "shade"),
                         min_pop_size = 1, color_period = 11, bar_height = 0.80,
                         base_colors = c("#e41a1c", "#377eb8", "#4daf4a",
                                         "#984ea3", "#ff7f00", "#ffff33")) {
  sort <- match.arg(sort)
  data <- .canonicalize(data, cols)
  core <- .prep_core(data, bar_height, arrange_by_frac = (sort == "frac"))

  prepped <- switch(
    sort,
    population = .prep_rainbow_colors(core, min_pop_size, color_period),
    frac       = .prep_rainbow_colors(core, min_pop_size, color_period),
    shade      = .prep_shade_colors(core, min_pop_size, base_colors)
  )

  .decanonicalize(prepped, cols)
}
