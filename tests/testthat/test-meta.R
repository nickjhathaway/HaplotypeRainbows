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

sample_meta <- function() {
  e <- environment()
  data("pfIsosHeomeV1_sampleMeta", package = "HaplotypeRainbows", envir = e)
  e$pfIsosHeomeV1_sampleMeta
}
target_meta <- function() {
  e <- environment()
  data("pfIsosHeomeV1_targetMeta", package = "HaplotypeRainbows", envir = e)
  e$pfIsosHeomeV1_targetMeta
}

test_that("set_sample_meta stores meta keyed by the sample column", {
  rb <- new_rb()$set_sample_meta(sample_meta(), match_col = "sample")
  sm <- rb$get_sample_meta()
  expect_true("sample" %in% names(sm))        # match_col renamed to the canonical key
  expect_true(all(c("country", "region") %in% names(sm)))
})

test_that("set_*_meta respects the cols argument", {
  rb <- new_rb()$set_sample_meta(sample_meta(), match_col = "sample",
                                 cols = c("country", "region"))
  expect_setequal(names(rb$get_sample_meta()), c("sample", "country", "region"))
})

test_that("set_meta errors on a missing match_col or column", {
  expect_error(new_rb()$set_sample_meta(sample_meta(), match_col = "nope"),
               "match_col")
  expect_error(
    new_rb()$set_sample_meta(sample_meta(), match_col = "sample", cols = "nope"),
    "not found"
  )
})

test_that("set_meta warns and NA-fills samples missing from metadata", {
  meta <- sample_meta()
  meta <- meta[meta$sample != "3D7", , drop = FALSE]   # drop one sample
  expect_warning(new_rb()$set_sample_meta(meta, match_col = "sample"),
                 "missing from the metadata")
})

test_that("add = TRUE merges columns; default replaces", {
  m <- sample_meta()
  rb <- new_rb()$set_sample_meta(m, match_col = "sample", cols = "country")
  rb$set_sample_meta(m, match_col = "sample", cols = "region", add = TRUE)
  expect_true(all(c("country", "region") %in% names(rb$get_sample_meta())))
  # default replace drops the earlier column
  rb$set_sample_meta(m, match_col = "sample", cols = "region")
  expect_false("country" %in% names(rb$get_sample_meta()))
})

test_that("sort_samples_by_meta reorders samples by metadata", {
  rb <- new_rb()$prep(sort = "population_rank")$
    set_sample_meta(sample_meta(), match_col = "sample")$
    sort_samples_by_meta("country")
  lv <- levels(rb$get_prepped()$sample)
  sm <- rb$get_sample_meta()
  countries <- sm$country[match(lv, sm$sample)]
  expect_false(is.unsorted(as.integer(factor(countries)), na.rm = TRUE))
})

test_that("sort_targets_by_meta and set_target_order reorder targets", {
  rb <- new_rb()$prep(sort = "population_rank")$
    set_target_meta(target_meta(), match_col = "target_name")$
    sort_targets_by_meta("class")
  expect_s3_class(rb$get_prepped()$target, "factor")

  first_two <- head(levels(rb$get_prepped()$target), 2)
  rb$set_target_order(rev(first_two))
  expect_identical(head(levels(rb$get_prepped()$target), 2), rev(first_two))
})

test_that("sort_by_meta errors without metadata set", {
  expect_error(new_rb()$prep(sort = "population_rank")$sort_samples_by_meta("country"),
               "No sample metadata")
})
