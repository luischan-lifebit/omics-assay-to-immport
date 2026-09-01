#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# salmon_to_immport_rnaseq.R  (Nextflow-ready version)
#
# Converts merged Salmon RNA-seq TPM matrices into ImmPort's RNA_SEQ_Results
# long format (schema version 3.37).
#
# Unlike the original script, paths and config are now passed as CLI args
# so this can run as a Nextflow process without editing the file.
#
# Usage:
#   salmon_to_immport_rnaseq.R <gene_tpm.tsv> <transcript_tpm.tsv> \
#       [--repository_name Ensembl] [--transcript_type mRNA] [--outdir .]
#
# Output (written to --outdir, default: current directory):
#   RNA_SEQ_Results_gene.tsv
#   RNA_SEQ_Results_transcript.tsv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(optparse)
})

# ---- CLI ARGS ---------------------------------------------------------------

option_list <- list(
  make_option("--repository_name", type = "character", default = "Ensembl",
              help = "Repository Name value [default: %default]"),
  make_option("--transcript_type", type = "character", default = "mRNA",
              help = "Transcript Type Reported value, applied to all rows [default: %default]"),
  make_option("--result_unit", type = "character", default = "TPM",
              help = "Result Unit Reported value [default: %default]"),
  make_option("--outdir", type = "character", default = ".",
              help = "Output directory [default: current dir]")
)

parser <- OptionParser(
  usage = "%prog <gene_tpm.tsv> <transcript_tpm.tsv> [options]",
  option_list = option_list
)
parsed <- parse_args(parser, positional_arguments = 2)
opts <- parsed$options
gene_counts_path <- parsed$args[1]
tx_counts_path   <- parsed$args[2]

dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- CONTROLLED VOCAB (from RNA_SEQ_Results.xlsx "lookup" sheet) ----------

valid_repository_names <- c(
  "ArrayExpress", "Broad Single Cell Portal", "dbGAP", "ENA", "Ensembl",
  "FlowRepository", "GenBank", "GEO", "GISAID", "IEDB", "ImmPort",
  "MassIVE", "MetaboLights", "Metabolomics Workbench", "MGnify",
  "NCBI Gene", "PRIDE", "SRA", "UniProt"
)
valid_transcript_types <- c("lincRNA", "mRNA", "snRNA")
valid_result_units     <- c("FPKM", "Gy", "Not Specified", "RPKM", "TPM")

stopifnot(opts$repository_name %in% valid_repository_names)
stopifnot(opts$transcript_type %in% valid_transcript_types)
stopifnot(opts$result_unit %in% valid_result_units)

# ---- SAMPLE MAP -------------------------------------------------------------
# Sample names are derived directly from the input file's column headers,
# passed through unchanged as Expsample IDs (still a placeholder until real
# digested ImmPort Expsample IDs exist -- see README open items).

derive_sample_map <- function(path, n_id_cols) {
  header <- names(read_tsv(path, n_max = 0, show_col_types = FALSE))
  samples <- header[(n_id_cols + 1):length(header)]
  setNames(samples, samples)
}

# ---- CORE CONVERSION FUNCTION ----------------------------------------------

melt_to_immport <- function(counts_path, id_col, sample_map) {
  df <- read_tsv(counts_path, show_col_types = FALSE)

  sample_cols <- names(sample_map)
  missing <- setdiff(sample_cols, colnames(df))
  if (length(missing) > 0) {
    stop("Sample columns not found in ", counts_path, ": ",
         paste(missing, collapse = ", "))
  }

  df %>%
    select(all_of(id_col), all_of(sample_cols)) %>%
    rename(`Reference Transcript ID` = !!id_col) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "Value Reported"
    ) %>%
    mutate(
      `Expsample ID`             = sample_map[sample],
      `Repository Name`          = opts$repository_name,
      `Transcript Type Reported` = opts$transcript_type,
      `Result Unit Reported`     = opts$result_unit,
      `Comments`                 = ""
    ) %>%
    select(
      `Expsample ID`,
      `Reference Transcript ID`,
      `Repository Name`,
      `Transcript Type Reported`,
      `Result Unit Reported`,
      `Value Reported`,
      `Comments`
    )
}

# ---- RUN ---------------------------------------------------------------------

gene_sample_map <- derive_sample_map(gene_counts_path, n_id_cols = 2)  # gene_id, gene_name
tx_sample_map   <- derive_sample_map(tx_counts_path, n_id_cols = 2)    # tx, gene_id

gene_results <- melt_to_immport(gene_counts_path, "gene_id", gene_sample_map)
gene_out <- file.path(opts$outdir, "RNA_SEQ_Results_gene.tsv")
write_tsv(gene_results, gene_out)
cat("Wrote", nrow(gene_results), "rows to", gene_out, "\n")

tx_results <- melt_to_immport(tx_counts_path, "tx", tx_sample_map)
tx_out <- file.path(opts$outdir, "RNA_SEQ_Results_transcript.tsv")
write_tsv(tx_results, tx_out)
cat("Wrote", nrow(tx_results), "rows to", tx_out, "\n")

# ---- SUMMARY -----------------------------------------------------------------

cat("\n--- CONFIRM BEFORE SUBMISSION ---\n")
cat("1. Repository Name = '", opts$repository_name,
    "' -- confirm correct for this data's ID system.\n", sep = "")
cat("2. Transcript Type Reported = '", opts$transcript_type,
    "' applied to ALL rows -- needs real biotype data to be fully accurate.\n", sep = "")
cat("3. Expsample ID = raw sample names from input file headers -- swap for\n",
    "   real digested ImmPort Expsample IDs once available.\n", sep = "")
