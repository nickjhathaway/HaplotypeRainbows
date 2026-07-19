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

test_that("rank_colors plots discretely and get_rank_colors maps invariant", {
  rb <- new_rb()$prep(sort = "population_rank", mark_invariant = TRUE)
  expect_s3_class(suppressWarnings(rb$plot(rank_colors = TRUE)), "ggplot")
  # n_ranks within the palette -> ranks map straight to the palette
  m <- rb$get_rank_colors(n_ranks = 5)
  expect_true("invariant" %in% names(m))
  expect_length(m, 6)
  expect_identical(unname(m["1"]), colorPalette_12[[1]])
})

test_that("rank_colors requires a rainbow prep, not shade", {
  expect_error(
    new_rb()$prep_shade()$plot(rank_colors = TRUE), "rank_colors requires"
  )
})

test_that("rank colours interpolate a ramp when ranks exceed the palette", {
  rb <- new_rb()$prep(sort = "population_rank")
  n <- max(rb$get_prepped()$popid)
  small <- colorPalette_08[1:3]
  # more ranks than the 3-colour palette -> ramp + warning (muffle the unrelated
  # ggplot "unknown aesthetics" warning so only the ramp warning is asserted)
  expect_warning(
    withCallingHandlers(
      rb$plot(rank_colors = TRUE, rank_palette = small),
      warning = function(w) {
        if (grepl("unknown aesthetic", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    "interpolating a colour ramp"
  )
  m <- rb$get_rank_colors(rank_palette = small)
  expect_length(m, n + 1)   # n ranks + invariant
})

test_that(".expand_rank_palette drops white/black and ramps", {
  pal <- c("#FF0000", "#00FF00", "#0000FF", "#FFFFFF", "#000000")
  out <- HaplotypeRainbows:::.expand_rank_palette(pal, 10)
  expect_length(out, 10)
  expect_false(any(toupper(out) %in% c("#FFFFFF", "#000000")))
  # when n <= palette length, returns the first n unchanged
  expect_identical(HaplotypeRainbows:::.expand_rank_palette(pal, 2), pal[1:2])
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
