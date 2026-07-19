# HaplotypeRainbows 2.0.1

## Bug fixes

* `prep()` / `prep_shade()` no longer error with "Names must be unique" when the input
  data contains a column whose name clashes with an internal column (e.g. a leftover
  SeekDeep `p_name` kept in the data after the target was remapped to another column).
  The engine now selects only the four mapped columns from the input and drops the rest.

## Internal changes

* The prepped table now uses a small, fixed, canonical column schema regardless of the
  input column names: `sample`, `target`, `hapid` and `rel_abund`, plus snake_case
  derived columns. (Previously it carried the user's original column names and the
  SeekDeep-style internal names.) `get_prepped()` therefore returns these canonical
  names; pass `get_prepped(original_names = TRUE)` to relabel the four key columns back
  to your input names.
* `rel_abund` now keeps the **raw** counts. The within-sample/target fraction (used for
  the bar geometry) is exposed as the new column `within_sample_freq`, with
  `total_abund` the per-(sample, target) denominator.

## New features

* `column_map()` returns the key mapping the canonical column names to the original
  input column names, so prepped data can be relabelled and rejoined to the source data.
* `get_prepped(original_names = TRUE)` relabels the canonical `sample` / `target` /
  `hapid` / `rel_abund` columns back to the input column names.
