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

test_that("plot_gap and legend column control are accepted", {
  rb <- prepped_rb()
  p <- suppressWarnings(rb$plot())
  expect_s3_class(
    suppressWarnings(rb$add_sample_metadata(
      p, cols = c("country", "region"), plot_gap = 2,
      legend_ncol = c(country = 4, region = 2)
    )),
    "ggplot"
  )
  expect_s3_class(
    suppressWarnings(rb$add_target_annotation(
      p, cols = "class", plot_gap = 1, legend_nrow = 1
    )),
    "ggplot"
  )
})

test_that("level_order controls legend/level order", {
  df <- data.frame(x = c("b", "a", "c", "a"))
  al <- HaplotypeRainbows:::.apply_level_order(df, "x", list(x = c("c", "a")))
  expect_identical(levels(al$x), c("c", "a", "b"))
  # unchanged column when no order supplied -> sorted
  al2 <- HaplotypeRainbows:::.apply_level_order(df, "x", NULL)
  expect_identical(levels(al2$x), c("a", "b", "c"))

  rb <- prepped_rb()
  p <- suppressWarnings(rb$plot())
  expect_s3_class(
    suppressWarnings(rb$add_sample_metadata(
      p, cols = "country", level_order = list(country = c("Uganda", "Kenya"))
    )),
    "ggplot"
  )
})

test_that("expanded rank palette makes adjacent colours distinct", {
  ramp <- HaplotypeRainbows:::.expand_rank_palette(colorPalette_12, 16)
  smooth <- grDevices::colorRampPalette(colorPalette_12)(16)
  adj_min <- function(cols) {
    rgb <- t(grDevices::col2rgb(cols))
    min(sqrt(rowSums((rgb[-1, ] - rgb[-nrow(rgb), ])^2)))
  }
  expect_gt(adj_min(ramp), adj_min(smooth))
})

test_that(".per_col resolves NULL/scalar/named/positional", {
  expect_equal(HaplotypeRainbows:::.per_col(NULL, c("a", "b")),
               list(a = NULL, b = NULL))
  expect_equal(HaplotypeRainbows:::.per_col(3, c("a", "b")),
               list(a = 3, b = 3))
  expect_equal(HaplotypeRainbows:::.per_col(c(a = 2), c("a", "b")),
               list(a = 2, b = NULL))
  expect_equal(HaplotypeRainbows:::.per_col(c(5, 6), c("a", "b")),
               list(a = 5, b = 6))
})

test_that("legend export: extract_legend, save_legend_pdf, drop_legends", {
  rb <- prepped_rb()
  p <- suppressWarnings(rb$add_sample_metadata(suppressWarnings(rb$plot()),
                                               cols = "region"))
  expect_s3_class(rb$extract_legend(p), "gtable")

  f <- tempfile(fileext = ".pdf")
  suppressWarnings(rb$save_legend_pdf(p, f, device = "pdf", width = 5, height = 6))
  expect_true(file.exists(f))
  expect_gt(file.info(f)$size, 0)

  nl <- rb$drop_legends(p)
  expect_s3_class(nl, "ggplot")
  expect_identical(nl$theme$legend.position, "none")
})
