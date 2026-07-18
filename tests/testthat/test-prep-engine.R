# Golden tests: the R6 engine must reproduce the output of the original standalone
# functions (snapshotted into _fixtures/ before the refactor). The example data uses
# the historical SeekDeep column names, so we construct the object with those names.

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

# Compare two prepped frames independent of row order and column order.
expect_prep_equal <- function(got, expected) {
  keys <- c("s_Sample", "p_name", "h_popUID")
  norm <- function(df) {
    df <- as.data.frame(dplyr::ungroup(df))
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
  rb <- new_rb()$prep(sort = "population")
  expect_prep_equal(rb$get_prepped(), fixture("prep_population.rds"))
})

test_that("frac prep matches the original prepForRainbowArrangedByFrac", {
  rb <- new_rb()$prep(sort = "frac")
  expect_prep_equal(rb$get_prepped(), fixture("prep_frac.rds"))
})

test_that("shade prep matches the original prepForRainbowShade", {
  rb <- new_rb()$prep(sort = "shade", min_pop_size = 1)
  expect_prep_equal(rb$get_prepped(), fixture("prep_shade.rds"))
})

test_that("clustering reorder matches the original (population)", {
  rb <- new_rb()$prep(sort = "population")$sort_by_clustering()
  expect_identical(levels(rb$get_prepped()$s_Sample),
                   fixture("resort_population_levels.rds"))
})

test_that("clustering reorder matches the original (major allele)", {
  rb <- new_rb()$prep(sort = "population")$sort_by_clustering(by_major_allele = TRUE)
  expect_identical(levels(rb$get_prepped()$s_Sample),
                   fixture("resort_major_levels.rds"))
})
