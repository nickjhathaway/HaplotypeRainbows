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
prepped_rb <- function() {
  new_rb()$prep(sort = "population_rank")$
    set_sample_meta(sample_meta(), match_col = "sample")$
    set_target_meta(target_meta(), match_col = "target_name")
}

test_that(".auto_meta_colors assigns one colour per distinct value", {
  df <- data.frame(a = c("x", "y", "z"), b = c("p", "p", "q"))
  cols <- HaplotypeRainbows:::.auto_meta_colors(df, c("a", "b"))
  expect_named(cols, c("a", "b"))
  expect_length(cols$a, 3)
  expect_length(cols$b, 2)
  expect_true(all(grepl("^#", cols$a)))
})

test_that(".auto_meta_colors interpolates beyond 15 levels", {
  df <- data.frame(a = as.character(1:20))
  cols <- HaplotypeRainbows:::.auto_meta_colors(df, "a")
  expect_length(cols$a, 20)
})

test_that("add_sample_metadata returns a ggplot", {
  rb <- prepped_rb()
  p <- suppressWarnings(rb$plot())
  p2 <- suppressWarnings(rb$add_sample_metadata(p, cols = c("country", "region")))
  expect_s3_class(p2, "ggplot")
})

test_that("add_target_annotation returns a ggplot", {
  rb <- prepped_rb()
  p <- suppressWarnings(rb$plot())
  p2 <- suppressWarnings(rb$add_target_annotation(p, cols = "class"))
  expect_s3_class(p2, "ggplot")
})

test_that("side and position are honoured", {
  rb <- prepped_rb()
  p <- suppressWarnings(rb$plot())
  expect_s3_class(
    suppressWarnings(rb$add_sample_metadata(p, cols = "country", side = "right")),
    "ggplot"
  )
  expect_s3_class(
    suppressWarnings(rb$add_target_annotation(p, cols = "class", position = "bottom")),
    "ggplot"
  )
})

test_that("annotation methods error without metadata / on bad columns", {
  rb <- new_rb()$prep(sort = "population_rank")
  p <- suppressWarnings(rb$plot())
  expect_error(rb$add_sample_metadata(p), "No sample metadata")
  expect_error(rb$add_target_annotation(p), "No target metadata")

  rb2 <- prepped_rb()
  expect_error(
    suppressWarnings(rb2$add_sample_metadata(suppressWarnings(rb2$plot()), cols = "nope")),
    "not found"
  )
})

test_that("colors override is applied", {
  rb <- prepped_rb()
  p <- suppressWarnings(rb$plot())
  override <- list(class = c(Diversity = "#111111", Drug = "#222222"))
  p2 <- suppressWarnings(
    rb$add_target_annotation(p, cols = "class", colors = override)
  )
  expect_s3_class(p2, "ggplot")
})
