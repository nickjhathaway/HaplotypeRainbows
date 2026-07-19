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

test_that("dims scales with target/sample counts and respects minimums", {
  rb <- new_rb()$prep(sort = "population_rank")
  d <- rb$dims(cell_width = 0.3, cell_height = 0.3)
  expect_true(d$width >= 6 && d$height >= 6)
  # 100 targets * 0.3 = 30 wide; 61 samples * 0.3 = 18.3 tall
  expect_equal(d$width, 30)
  expect_equal(d$height, 18.3)
  expect_equal(rb$dims(extra_width = 5)$width, 35)
})

test_that("save_pdf writes a file (base pdf device)", {
  rb <- new_rb()$prep(sort = "population_rank")
  p <- suppressWarnings(rb$plot())
  f <- tempfile(fileext = ".pdf")
  suppressWarnings(rb$save_pdf(p, f, width = 8, height = 6, device = "pdf"))
  expect_true(file.exists(f))
  expect_gt(file.info(f)$size, 0)
})
