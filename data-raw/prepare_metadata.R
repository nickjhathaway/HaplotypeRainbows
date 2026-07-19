# Prepare the example metadata datasets that accompany pfIsosHeomeV1.
# Source data (not shipped): the full biosample metadata table and the target insert
# table. We subset to the samples/targets present in the example data and drop columns
# that are constant across the subset.
#
# Run with: Rscript data-raw/prepare_metadata.R

library(dplyr)
library(readr)

load("data/pfisolateExample.rda")
samps <- unique(as.character(pfIsosHeomeV1$s_Sample))
tars  <- unique(as.character(pfIsosHeomeV1$p_name))

sample_src <- "/Users/nhathaway/Documents/tank/data/plasmodium/falciparum/pfdata/metadata/metaByBioSample.tab.txt"
target_src <- "/Users/nhathaway/Downloads/heome1_Pf3D7_inserts_out.tsv"

# --- sample metadata: subset to example samples, drop constant columns (keep `sample`)
pfIsosHeomeV1_sampleMeta <- read_tsv(sample_src, show_col_types = FALSE) %>%
  filter(sample %in% samps) %>%
  distinct(sample, .keep_all = TRUE)
constant_cols <- names(pfIsosHeomeV1_sampleMeta)[
  vapply(pfIsosHeomeV1_sampleMeta, function(x) dplyr::n_distinct(x) <= 1, logical(1))
]
constant_cols <- setdiff(constant_cols, "sample")
pfIsosHeomeV1_sampleMeta <- pfIsosHeomeV1_sampleMeta %>%
  select(-dplyr::all_of(constant_cols)) %>%
  relocate(sample) %>%
  as.data.frame()

# --- target metadata: subset to example targets
pfIsosHeomeV1_targetMeta <- read_tsv(target_src, show_col_types = FALSE) %>%
  filter(target_name %in% tars) %>%
  distinct(target_name, .keep_all = TRUE) %>%
  as.data.frame()

message("sample meta: ", nrow(pfIsosHeomeV1_sampleMeta), " rows, cols: ",
        paste(names(pfIsosHeomeV1_sampleMeta), collapse = ", "))
message("dropped constant sample cols: ", paste(constant_cols, collapse = ", "))
message("target meta: ", nrow(pfIsosHeomeV1_targetMeta), " rows, cols: ",
        paste(names(pfIsosHeomeV1_targetMeta), collapse = ", "))

save(pfIsosHeomeV1_sampleMeta, file = "data/pfIsosHeomeV1_sampleMeta.rda",
     compress = "xz")
save(pfIsosHeomeV1_targetMeta, file = "data/pfIsosHeomeV1_targetMeta.rda",
     compress = "xz")
