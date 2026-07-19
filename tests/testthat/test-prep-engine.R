# Golden tests: the R6 engine must reproduce the geometry / colours / ordering of the
# original standalone functions (snapshotted into _fixtures/ before the refactor). The
# v2.0.1 engine renames its columns to a canonical snake_case schema (sample, target,
# hapid, rel_abund + derived) and keeps rel_abund RAW (adding within_sample_freq /
# total_abund), so we map the new output back onto the fixtures' original schema before
# comparing. The example data uses the historical SeekDeep column names, so we construct
# the object with those names.

fixture <- function(name) {
  readRDS(testthat::test_path("_fixtures", name))
}

new_rb <- function() {
  data("pfisolateExample", package = "HaplotypeRainbows", envir = environment())
  HaplotypeRainbow$new(
    pfIsosHeomeV1,
    sample_col    = "s_Sample",
    target_col    = "p_name",
    popuid_col    = "h_popUID",
    rel_abund_col = "c_AveragedFrac"
  )
}

# Canonical/snake_case column names produced by the v2.0.1 engine -> the original names
# stored in the fixtures. Columns kept with the same name in both schemas (samp_n, h_id,
# h_id_freq, h_id_freq_mod, p_hue, n, total, freq, h_color, h_color_mod) are omitted.
.new_to_old <- c(
  sample = "s_Sample", target = "p_name", hapid = "h_popUID",
  within_sample_freq = "c_AveragedFrac", total_abund = "totalAbund",
  n_targets = "targetNumber", sample_coi = "s_COI",
  within_sample_freq_mod = "relAbundCol_mod",
  freq_cumsum = "fracCumSum", freq_mod_cumsum = "fracModCumSum",
  fake_freq = "fakeFrac", fake_freq_mod = "fakeFracMod",
  fake_freq_cumsum = "fakeFracCumSum", fake_freq_mod_cumsum = "fakeFracModCumSum",
  pop_id = "popid", max_pop_id = "maxPopid", pop_id_frac = "popidFrac",
  hue_mod = "hueMod", pop_id_perc = "popidPerc",
  pop_id_frac_reg_color = "popidFracRegColor", pop_id_perc_log = "popidPercLog",
  pop_id_frac_log_color = "popidFracLogColor",
  p_uniq_haps = "p_uniqHaps", cum_freq = "cumFreq", mod_cum_freq = "modCumFreq",
  p_color_id = "p_color_ID", p_base_color = "p_baseColor",
  h_color_by_freq = "h_color_byFreq", h_color_by_freq_mod = "h_color_byFreq_mod"
)

# The shade RGBA hex columns are a deterministic function of numeric columns we already
# compare (cumFreq, modCumFreq, h_id_freq, h_id_freq_mod, p_baseColor), but scales::alpha()
# quantizes the numeric alpha to an 8-bit byte and rounds differently across platforms at
# boundaries (e.g. macOS "#FF7F00B3" vs Linux "#FF7F00B2"). We drop these from the
# byte-exact comparison and verify the numeric drivers instead. Named in the fixtures'
# (original) schema, since we compare after mapping back to it.
.volatile_color_cols <- c("h_color", "h_color_mod",
                          "h_color_byFreq", "h_color_byFreq_mod")

# Map the new engine's output back onto the fixtures' original schema. The new raw
# `rel_abund` has no equivalent in the old fixtures (there c_AveragedFrac WAS the
# fraction), so it is dropped; within_sample_freq maps onto the old c_AveragedFrac.
to_old_schema <- function(df) {
  df <- as.data.frame(dplyr::ungroup(df))
  df[["rel_abund"]] <- NULL
  nm <- names(df)
  hit <- nm %in% names(.new_to_old)
  names(df)[hit] <- unname(.new_to_old[nm[hit]])
  df
}

# Compare two prepped frames independent of row order and column order.
expect_prep_equal <- function(got, expected) {
  got <- to_old_schema(got)
  keys <- c("s_Sample", "p_name", "h_popUID")
  norm <- function(df) {
    df <- as.data.frame(dplyr::ungroup(df))
    df <- df[, !names(df) %in% .volatile_color_cols, drop = FALSE]
    df <- df[do.call(order, lapply(df[keys], as.character)), , drop = FALSE]
    df <- df[, order(names(df)), drop = FALSE]
    # compare factor columns by value; their levels are asserted separately
    factor_cols <- vapply(df, is.factor, logical(1))
    df[factor_cols] <- lapply(df[factor_cols], as.character)
    rownames(df) <- NULL
    df
  }
  # factor levels (drive plotting order) must match too
  testthat::expect_identical(levels(got$s_Sample), levels(expected$s_Sample))
  testthat::expect_identical(levels(got$p_name), levels(expected$p_name))
  testthat::expect_equal(norm(got), norm(expected))
}

test_that("population prep matches the original prepForRainbow", {
  rb <- new_rb()$prep(sort = "population_rank")
  expect_prep_equal(rb$get_prepped(), fixture("prep_population.rds"))
})

test_that("frac prep matches the original prepForRainbowArrangedByFrac", {
  rb <- new_rb()$prep(sort = "within_sample_freq")
  expect_prep_equal(rb$get_prepped(), fixture("prep_frac.rds"))
})

test_that("shade prep matches the original prepForRainbowShade", {
  rb <- new_rb()$prep_shade(min_pop_size = 1)
  expect_prep_equal(rb$get_prepped(), fixture("prep_shade.rds"))
})

test_that("clustering reorder matches the original (population)", {
  rb <- new_rb()$prep(sort = "population_rank")$sort_by_clustering()
  expect_identical(levels(rb$get_prepped()$sample),
                   fixture("resort_population_levels.rds"))
})

test_that("clustering reorder matches the original (major allele)", {
  rb <- new_rb()$prep(sort = "population_rank")$sort_by_clustering(by_major_allele = TRUE)
  expect_identical(levels(rb$get_prepped()$sample),
                   fixture("resort_major_levels.rds"))
})

test_that("prepped table uses the canonical schema; rel_abund stays raw", {
  p <- dplyr::ungroup(new_rb()$prep(sort = "population_rank")$get_prepped())
  # canonical names present, original SeekDeep names gone
  expect_true(all(c("sample", "target", "hapid", "rel_abund",
                    "within_sample_freq", "total_abund") %in% names(p)))
  expect_false(any(c("s_Sample", "p_name", "h_popUID", "c_AveragedFrac") %in% names(p)))
  # within_sample_freq is rel_abund normalised within each (sample, target)
  expect_equal(p$within_sample_freq, p$rel_abund / p$total_abund)
  chk <- p %>%
    dplyr::group_by(sample, target) %>%
    dplyr::summarise(s = sum(rel_abund), t = dplyr::first(total_abund),
                     f = sum(within_sample_freq), .groups = "drop")
  expect_equal(chk$s, chk$t)        # total_abund is the per-cell raw sum
  expect_equal(chk$f, rep(1, nrow(chk)))  # fractions sum to 1 per cell
})

test_that("prep() drops non-mapped columns, avoiding canonical-name collisions", {
  # Mirrors the real bug: the target is remapped to a new column while a stray column
  # named like a canonical internal name (here literally "target"/"sample") lingers.
  df <- data.frame(
    s_Sample       = c("A", "A", "B"),
    my_target      = c("t1", "t1", "t1"),
    h_popUID       = c("h1", "h2", "h1"),
    c_AveragedFrac = c(5, 3, 8),
    target         = c("STRAY", "STRAY", "STRAY"),  # clashes with canonical name
    sample         = c("X", "X", "Y"),              # clashes with canonical name
    stringsAsFactors = FALSE
  )
  rb <- HaplotypeRainbow$new(df, sample_col = "s_Sample", target_col = "my_target",
                             popuid_col = "h_popUID", rel_abund_col = "c_AveragedFrac")
  expect_error(rb$prep(), NA)                 # no "Names must be unique" collision
  p <- rb$get_prepped()
  expect_setequal(unique(as.character(p$target)), "t1")   # the real target, not STRAY
})

test_that("column_map / get_prepped(original_names) reconnect to the input names", {
  rb <- new_rb()$prep(sort = "population_rank")
  cm <- rb$column_map()
  expect_identical(unname(cm[c("sample", "target", "hapid", "rel_abund")]),
                   c("s_Sample", "p_name", "h_popUID", "c_AveragedFrac"))
  orig <- rb$get_prepped(original_names = TRUE)
  expect_true(all(c("s_Sample", "p_name", "h_popUID", "c_AveragedFrac") %in% names(orig)))
  expect_false(any(c("sample", "target", "hapid", "rel_abund") %in% names(orig)))
})
