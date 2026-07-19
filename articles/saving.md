# Saving, sizing & interactivity

Rainbow plots get wide and tall fast, so the class computes figure sizes
for you and writes PDFs sized to the data. This article covers
auto-sizing, saving, exporting the legend on its own, and interactive
tooltips.

``` r

library(HaplotypeRainbows)
library(ggplot2)

data("pfisolateExample")
rb <- HaplotypeRainbow$new(
  pfIsosHeomeV1,
  sample_col    = "s_Sample",
  target_col    = "p_name",
  popuid_col    = "h_popUID",
  rel_abund_col = "c_AveragedFrac"
)
rb$prep(sort = "population_rank")
```

## Figure size from the data

`dims()` scales width with the number of targets and height with the
number of samples (counting any inter-cluster spacer rows), respecting
minimums. Use it to size a Quarto / R Markdown chunk, or let
`save_pdf()` use it automatically.

``` r

rb$dims()                                  # list(width, height) in inches
#> $width
#> [1] 30
#> 
#> $height
#> [1] 18.3
rb$dims(cell_width = 0.4, cell_height = 0.4, extra_width = 4)
#> $width
#> [1] 44
#> 
#> $height
#> [1] 24.4
```

## Saving to PDF

`save_pdf()` sizes the figure from `dims()` (unless you pass
`width`/`height`) and uses `cairo_pdf` by default, which keeps hyphens
as hyphens (falling back to
[`pdf()`](https://rdrr.io/r/grDevices/pdf.html) when the platform lacks
cairo). The background is transparent by default.

``` r

p <- rb$plot()
rb$save_pdf(p, "rainbow.pdf")                  # auto-sized, cairo device
rb$save_pdf(p, "rainbow.pdf", device = "pdf")  # base pdf device (e.g. Windows)
rb$save_pdf(p, "rainbow.pdf", bg = "white")    # opaque background

# reserve room for the legend in the page size too (off by default):
rb$save_pdf(p, "rainbow.pdf", size_include_legend = TRUE)
```

## Exporting the legend separately

Legends on wide metadata can be large. To place the plot and its legend
independently in a figure, export them separately: `save_legend_pdf()`
writes **just** the legend (auto-sized so it isn’t clipped),
`drop_legends()` returns the plot without it, and `extract_legend()`
returns the legend grob for composing with patchwork / cowplot / grid.

``` r

rb$save_legend_pdf(p, "legend.pdf")             # the legend only
rb$save_pdf(rb$drop_legends(p), "plot.pdf")     # the plot without its legend
legend_grob <- rb$extract_legend(p)             # or compose it yourself
```

## Interactive tooltips

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) carries the
sample, target, haplotype and abundance values as extra aesthetics, so
[`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html)
surfaces them on hover:

``` r

p <- rb$plot()
plotly::ggplotly(p)
```
