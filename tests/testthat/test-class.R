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

test_that("constructor errors on missing columns", {
  expect_error(
    HaplotypeRainbow$new(data.frame(a = 1)),
    "Column\\(s\\) not found"
  )
})

test_that("operations before prep() error clearly", {
  expect_error(new_rb()$plot(), "Call \\$prep")
  expect_error(new_rb()$sort_by_clustering(), "Call \\$prep")
})

test_that("prep/sort return self (chainable) and plot returns a ggplot", {
  rb <- new_rb()
  expect_identical(rb$prep(sort = "population_rank"), rb)
  expect_identical(rb$sort_alphabetical(), rb)
  expect_s3_class(suppressWarnings(rb$plot()), "ggplot")
})

test_that("rainbow style rejects shade-prepped data and vice versa", {
  expect_error(
    new_rb()$prep_shade()$plot(style = "rainbow"), "prep\\(\\)"
  )
  expect_error(
    new_rb()$prep(sort = "population_rank")$plot(style = "shade"), "prep_shade"
  )
})

test_that("x/y axis labels can be toggled independently", {
  rb <- new_rb()$prep(sort = "population_rank")
  expect_s3_class(suppressWarnings(
    rb$plot(x_axis_labels = FALSE, y_axis_labels = FALSE)
  ), "ggplot")
  expect_s3_class(suppressWarnings(
    rb$plot(x_axis_labels = TRUE, y_axis_labels = FALSE)
  ), "ggplot")
})

test_that("shade prep plots with identity fill by default", {
  p <- suppressWarnings(new_rb()$prep_shade(min_pop_size = 1)$plot())
  expect_s3_class(p, "ggplot")
})

test_that("set_sample_order puts requested samples first", {
  rb <- new_rb()$prep(sort = "population_rank")
  wanted <- c("3D7", "HB3")
  rb$set_sample_order(wanted)
  expect_identical(head(levels(rb$get_prepped()$s_Sample), 2), wanted)
})

test_that("haplotype_rainbow() wrapper builds the same class", {
  rb <- haplotype_rainbow(
    data.frame(library_sample_name = "s", target_name = "t",
               seq = "a", reads = 1)
  )
  expect_s3_class(rb, "HaplotypeRainbow")
})
