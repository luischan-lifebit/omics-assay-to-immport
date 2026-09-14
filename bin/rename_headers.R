#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# rename_headers.R
#
# Generic column-header renamer. Takes any TSV and a mapping CSV (old_name,
# new_name), and writes a copy with headers renamed accordingly. Used to
# produce ImmPort's canonical submission format (Title Case headers) from
# the lowercase_with_underscores format used internally / for OMOP ETL --
# but the script itself knows nothing about RNA-seq, ImmPort, or OMOP. It's
# generic on purpose, so the same script works for proteomics, demographics,
# or any future assay type: just supply a different mapping CSV.
#
# Usage:
#   rename_headers.R <input.tsv> <header_map.csv> <output.tsv>
#
# header_map.csv must have columns: old_name,new_name
# Any column in <input.tsv> not listed in the map is passed through
# unchanged (its name is kept as-is), rather than dropped or erroring --
# this keeps the script safe to reuse even if a mapping file is incomplete.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: rename_headers.R <input.tsv> <header_map.csv> <output.tsv>")
}
input_path  <- args[1]
map_path    <- args[2]
output_path <- args[3]

if (!file.exists(input_path)) stop("Input file not found: ", input_path)
if (!file.exists(map_path))   stop("Header map file not found: ", map_path)

header_map_df <- read_csv(map_path, show_col_types = FALSE)
required_cols <- c("old_name", "new_name")
missing_cols  <- setdiff(required_cols, colnames(header_map_df))
if (length(missing_cols) > 0) {
  stop("Header map is missing required column(s): ",
       paste(missing_cols, collapse = ", "))
}

header_map <- setNames(header_map_df$new_name, header_map_df$old_name)

df <- read_tsv(input_path, show_col_types = FALSE)

new_names <- header_map[colnames(df)]
# columns not present in the map keep their original name, unchanged
new_names[is.na(new_names)] <- colnames(df)[is.na(new_names)]
colnames(df) <- unname(new_names)

write_tsv(df, output_path)
cat("Wrote", nrow(df), "rows to", output_path, "\n")
