#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# salmon_to_immport_rnaseq.R  (Nextflow-ready, with optional linkage support)
#
# Converts merged Salmon RNA-seq TPM matrices into ImmPort's RNA_SEQ_Results
# long format (schema version 3.37).
#
# Usage:
#   salmon_to_immport_rnaseq.R <gene_tpm.tsv> <transcript_tpm.tsv> \
#       [--repository_name Ensembl] [--transcript_type mRNA] \
#       [--result_unit TPM] [--outdir .] [--linkage_file linkage.csv]
#
# --linkage_file is OPTIONAL. When provided, it must be a CSV with columns
# `source_person_id,sample_id` (per Lifebit's lifebit_omics_linkage contract:
# https://lifebit.atlassian.net/wiki/spaces/DEL/pages/2595454986)
# where sample_id matches this run's Expsample ID (our sample column headers)
# and source_person_id is the real participant identifier (must match
# person.person_source_value in the target OMOP schema downstream).
#
# When --linkage_file is supplied, a `Participant ID` column is added to the
# output, joined from sample_id -> source_person_id. Multiple samples can map
# to the same participant (e.g. baseline + follow-up), matching the real
# many-samples-to-one-participant pattern documented in the linkage contract.
#
# When --linkage_file is omitted, output is unchanged from before -- no
# Participant ID column is added.
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
})

# ---- CLI ARGS (base R only -- no optparse, no runtime install needed) -----

raw_args <- commandArgs(trailingOnly = TRUE)

is_flag <- grepl("^--", raw_args)
first_flag_idx <- which(is_flag)[1]

if (is.na(first_flag_idx)) {
  positional <- raw_args
  flag_args  <- character(0)
} else {
  positional <- raw_args[seq_len(first_flag_idx - 1)]
  flag_args  <- raw_args[first_flag_idx:length(raw_args)]
}

if (length(positional) < 2) {
  stop("Usage: salmon_to_immport_rnaseq.R <gene_tpm.tsv> <transcript_tpm.tsv> ",
       "[--repository_name X] [--transcript_type X] [--result_unit X] ",
       "[--outdir X] [--linkage_file X]")
}
gene_counts_path <- positional[1]
tx_counts_path   <- positional[2]

parse_flags <- function(flag_args, defaults) {
  opts <- defaults
  i <- 1
  while (i <= length(flag_args)) {
    key <- sub("^--", "", flag_args[i])
    if (key %in% names(defaults) && i + 1 <= length(flag_args)) {
      opts[[key]] <- flag_args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  opts
}

opts <- parse_flags(flag_args, list(
  repository_name = "Ensembl",
  transcript_type  = "mRNA",
  result_unit      = "TPM",
  outdir           = ".",
  linkage_file     = NA_character_
))

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

# ---- OPTIONAL LINKAGE FILE --------------------------------------------------
# Per lifebit_omics_linkage contract: source_person_id,sample_id
# sample_id here must match our Expsample ID (sample column headers).

has_linkage <- !is.na(opts$linkage_file) && nzchar(opts$linkage_file)

linkage_map <- NULL
if (has_linkage) {
  if (!file.exists(opts$linkage_file)) {
    stop("--linkage_file provided but not found: ", opts$linkage_file)
  }
  linkage_df <- read_csv(opts$linkage_file, show_col_types = FALSE)

  required_cols <- c("source_person_id", "sample_id")
  missing_cols <- setdiff(required_cols, colnames(linkage_df))
  if (length(missing_cols) > 0) {
    stop("--linkage_file is missing required column(s): ",
         paste(missing_cols, collapse = ", "),
         ". Expected columns: source_person_id, sample_id")
  }

  # named vector: sample_id -> source_person_id (one participant can repeat
  # across multiple sample_id keys -- that's expected and fine)
  linkage_map <- setNames(linkage_df$source_person_id, linkage_df$sample_id)

  cat("Loaded linkage file:", opts$linkage_file,
      "--", nrow(linkage_df), "sample-to-participant mappings\n")
}

# ---- SAMPLE MAP -------------------------------------------------------------

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

  result <- df %>%
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
    )

  if (has_linkage) {
    result <- result %>%
      mutate(`Participant ID` = linkage_map[`Expsample ID`])

    unmatched <- result %>%
      filter(is.na(`Participant ID`)) %>%
      pull(`Expsample ID`) %>%
      unique()
    if (length(unmatched) > 0) {
      warning("No linkage entry found for ", length(unmatched),
              " Expsample ID(s): ", paste(unmatched, collapse = ", "),
              " -- Participant ID left blank for these rows.")
    }

    result <- result %>%
      select(
        `Participant ID`,
        `Expsample ID`,
        `Reference Transcript ID`,
        `Repository Name`,
        `Transcript Type Reported`,
        `Result Unit Reported`,
        `Value Reported`,
        `Comments`
      )
  } else {
    result <- result %>%
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

  result
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
if (has_linkage) {
  cat("3. Participant ID populated from linkage file '", opts$linkage_file,
      "'. Confirm source_person_id values match person.person_source_value\n",
      "   in the target OMOP schema before downstream ingestion.\n", sep = "")
} else {
  cat("3. No --linkage_file provided -- output has no Participant ID column.\n",
      "   Provide one (source_person_id,sample_id CSV) if this output feeds\n",
      "   into OMOP ingestion requiring participant linkage.\n", sep = "")
}
