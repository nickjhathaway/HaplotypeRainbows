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

test_that("get_hclust / get_dendrogram error before clustering, work after", {
  rb <- new_rb()$prep(sort = "population_rank")
  expect_error(rb$get_hclust(), "No clustering")
  rb$sort_by_clustering()
  expect_s3_class(rb$get_hclust(), "hclust")
  expect_s3_class(rb$get_dendrogram(), "dendrogram")
})

test_that("cluster_groups cuts into k groups, relabelled in dendrogram order", {
  rb <- new_rb()$prep(sort = "population_rank")$sort_by_clustering()
  cg <- rb$cluster_groups(k = 4)
  expect_true(all(c("s_Sample", "cluster") %in% names(cg)))
  expect_equal(nlevels(cg$cluster), 4)
  # first sample belongs to cluster 1 (labels are relabelled by appearance)
  expect_identical(as.character(cg$cluster[1]), "1")
  expect_error(rb$cluster_groups(), "Supply either k or h")
})

sample_meta <- function() {
  e <- environment()
  data("pfIsosHeomeV1_sampleMeta", package = "HaplotypeRainbows", envir = e)
  e$pfIsosHeomeV1_sampleMeta
}

test_that("export_groups_pdf combines per-cluster pages into one PDF", {
  skip_if_not_installed("qpdf")
  rb <- new_rb()$prep(sort = "population_rank")$sort_by_clustering()
  f <- tempfile(fileext = ".pdf")
  suppressWarnings(rb$export_groups_pdf(
    f, plot_fun = function(sub) suppressWarnings(sub$plot(x_axis_labels = FALSE)),
    by = "cluster", k = 4, device = "pdf"
  ))
  expect_true(file.exists(f))
  expect_equal(qpdf::pdf_length(f), 4L)
})

test_that("export_groups_pdf by meta can write separate files", {
  rb <- new_rb()$prep(sort = "population_rank")$
    set_sample_meta(sample_meta(), match_col = "sample")
  outs <- suppressWarnings(rb$export_groups_pdf(
    tempfile(fileext = ".pdf"),
    plot_fun = function(sub) suppressWarnings(sub$plot(x_axis_labels = FALSE)),
    by = "meta", meta_col = "secondaryRegion", combine = FALSE, device = "pdf"
  ))
  expect_gt(length(outs), 1)
  expect_true(all(file.exists(outs)))
})

test_that("export_groups_pdf validates grouping inputs", {
  rb <- new_rb()$prep(sort = "population_rank")
  expect_error(
    rb$export_groups_pdf(tempfile(), function(s) s$plot(),
                         by = "meta", meta_col = "x"),
    "No sample metadata"
  )
})

test_that("add_cluster_gaps inserts spacer rows between clusters", {
  rb <- new_rb()$prep(sort = "population_rank")$sort_by_clustering()
  n_real <- dplyr::n_distinct(rb$get_prepped()$s_Sample)
  rb$add_cluster_gaps(k = 6, gap = 2)
  # 6 clusters -> 5 gaps x 2 rows = 10 spacer levels added
  expect_equal(nlevels(rb$get_prepped()$s_Sample), n_real + 10)
  # plot renders and drops spacer labels
  expect_s3_class(suppressWarnings(rb$plot()), "ggplot")
  # gap = 0 removes the spacers
  rb$add_cluster_gaps(k = 6, gap = 0)
  expect_equal(nlevels(rb$get_prepped()$s_Sample), n_real)
})

test_that("add_cluster_gaps requires a prior clustering", {
  expect_error(
    new_rb()$prep(sort = "population_rank")$add_cluster_gaps(k = 3),
    "No clustering"
  )
})

test_that("add_cluster_bands overlays one band per group", {
  rb <- new_rb()$prep(sort = "population_rank")$sort_by_clustering()
  p <- suppressWarnings(rb$plot())
  expect_s3_class(suppressWarnings(rb$add_cluster_bands(p, k = 5)), "ggplot")
  expect_s3_class(
    suppressWarnings(rb$add_cluster_bands(p, k = 5, expand = 2, border = "black")),
    "ggplot"
  )
  # needs a k/h, and needs clustering first
  expect_error(rb$add_cluster_bands(p), "Supply either k or h")
  expect_error(
    new_rb()$prep(sort = "population_rank")$add_cluster_bands(p, k = 3),
    "No clustering"
  )
})

test_that("clustering honours dist/linkage and a target subset", {
  rb <- new_rb()$prep(sort = "population_rank")
  tars <- head(unique(as.character(rb$get_prepped()$p_name)), 10)
  rb$sort_by_clustering(targets = tars, dist_method = "manhattan",
                        hclust_method = "complete")
  expect_s3_class(rb$get_hclust(), "hclust")
  expect_identical(rb$get_hclust()$method, "complete")

  expect_error(
    new_rb()$prep(sort = "population_rank")$sort_by_clustering(targets = "nope"),
    "None of the supplied"
  )
})
