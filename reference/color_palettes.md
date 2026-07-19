# Colour-blind-friendly categorical palettes

Three categorical palettes (8, 12 and 15 colours) drawn from
colour-blind-friendly colour sets, suitable for colouring sample/target
metadata bands or for use as the `colors` ramp in
[HaplotypeRainbow](https://nickjhathaway.github.io/HaplotypeRainbows/reference/HaplotypeRainbow.md)'s
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Usage

``` r
colorPalette_08

colorPalette_12

colorPalette_15
```

## Format

Character vectors of hex colour strings:

- colorPalette_08:

  8 colours.

- colorPalette_12:

  12 colours.

- colorPalette_15:

  15 colours.

## Examples

``` r
scales::show_col(colorPalette_12)
```
