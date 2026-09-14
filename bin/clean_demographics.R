#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# clean_demographics.R
#
# Cleans a user-provided demographic file's column headers so it's directly
# compatible with Data Factory OMOP ETL. Does NOT generate or infer missing
# demographic data -- only reformats what the user provides.
#
# Required fields (per the OMOP person table's actual Required validation --
# confirmed directly against the Edit Person Table screen, not assumed):
#   participant_id, gender, year_of_birth, race, ethnicity
#
# NOTE: year_of_birth is the required date field, not date_of_birth /
# birth_datetime -- birth_datetime is supported by OMOP but not Required.
# date_of_birth is passed through if present, but its absence does not
# fail the workflow.
#
# Column matching is exact-label only (case/whitespace/underscore
# normalised), same convention as the RNA-seq script -- e.g. "Participant ID"
# and "participant id" both match participant_id, but a synonym like "Sex"
# or "DOB" does not auto-match gender / year_of_birth. Extend
# KNOWN_HEADER_ALIASES below if synonym matching is confirmed as needed.
#
# Usage:
#   clean_demographics.R <demographics_file> \
#       [--outdir .] [--linkage_file linkage.csv]
#
# --linkage_file is OPTIONAL. When provided (same lifebit_omics_linkage
# contract as the RNA-seq script: source_person_id,sample_id), this script
# checks whether the demographic file's participant_id values overlap with
# the linkage file's source_person_id values, and WARNS (does not fail) on
# any mismatch in either direction -- per AC4, this is currently implemented
# as a validation warning, pending confirmation from Sangram on whether it
# should hard-fail instead.
#
# Output (written to --outdir, default: current directory):
#   cleaned_demographics.tsv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
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

if (length(positional) < 1) {
  stop("Usage: clean_demographics.R <demographics_file> [--outdir X] [--linkage_file X]")
}
demographics_path <- positional[1]

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
  outdir       = ".",
  linkage_file = NA_character_
))

dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(demographics_path)) {
  stop("Demographics file not found: ", demographics_path)
}

# ---- REQUIRED / OPTIONAL FIELDS --------------------------------------------
# Confirmed against the OMOP person table's Required validation directly --
# year_of_birth is required, birth_datetime is not.

REQUIRED_FIELDS <- c("participant_id", "gender", "year_of_birth", "race", "ethnicity")
OPTIONAL_FIELDS <- c("date_of_birth")

# ---- HEADER NORMALISATION --------------------------------------------------
# Exact-label matching, case/whitespace/underscore-insensitive.
# "Participant ID", "participant id", "PARTICIPANT_ID" all normalise to the
# same key and match participant_id. Synonyms (Sex, DOB) do NOT match --
# see note at top of file.

normalise_header <- function(x) {
  x %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all("[\\s_]+", "_")
}

# ---- READ INPUT -------------------------------------------------------------

df <- read_delim(demographics_path, delim = NULL, show_col_types = FALSE)
original_headers <- colnames(df)
normalised_headers <- normalise_header(original_headers)

all_target_fields <- c(REQUIRED_FIELDS, OPTIONAL_FIELDS)
header_lookup <- setNames(original_headers, normalised_headers)

# ---- VALIDATE REQUIRED FIELDS (AC3: fail with clear warning if missing) ---

missing_required <- REQUIRED_FIELDS[!REQUIRED_FIELDS %in% normalised_headers]
if (length(missing_required) > 0) {
  stop(
    "Demographics file is missing required column(s): ",
    paste(missing_required, collapse = ", "),
    ". Required fields are: ", paste(REQUIRED_FIELDS, collapse = ", "),
    ". Found columns: ", paste(original_headers, collapse = ", ")
  )
}

# ---- BUILD CLEANED OUTPUT --------------------------------------------------

present_target_fields <- all_target_fields[all_target_fields %in% normalised_headers]
source_cols <- header_lookup[present_target_fields]

cleaned <- df %>%
  select(all_of(unname(source_cols))) %>%
  setNames(present_target_fields)

# ---- OPTIONAL LINKAGE FILE CONSISTENCY CHECK (AC4) -------------------------
# Warns, does not fail -- pending confirmation on whether this should be a
# hard failure instead.

has_linkage <- !is.na(opts$linkage_file) && nzchar(opts$linkage_file)

if (has_linkage) {
  if (!file.exists(opts$linkage_file)) {
    stop("--linkage_file provided but not found: ", opts$linkage_file)
  }
  linkage_df <- read_csv(opts$linkage_file, show_col_types = FALSE)

  required_linkage_cols <- c("source_person_id", "sample_id")
  missing_linkage_cols <- setdiff(required_linkage_cols, colnames(linkage_df))
  if (length(missing_linkage_cols) > 0) {
    stop("--linkage_file is missing required column(s): ",
         paste(missing_linkage_cols, collapse = ", "))
  }

  demographics_ids <- unique(cleaned$participant_id)
  linkage_ids       <- unique(linkage_df$source_person_id)

  in_demographics_not_linkage <- setdiff(demographics_ids, linkage_ids)
  in_linkage_not_demographics <- setdiff(linkage_ids, demographics_ids)

  if (length(in_demographics_not_linkage) > 0) {
    warning(
      length(in_demographics_not_linkage),
      " participant_id(s) in the demographics file have no matching ",
      "source_person_id in the linkage file: ",
      paste(in_demographics_not_linkage, collapse = ", ")
    )
  }
  if (length(in_linkage_not_demographics) > 0) {
    warning(
      length(in_linkage_not_demographics),
      " source_person_id(s) in the linkage file have no matching ",
      "participant_id in the demographics file -- those participants' ",
      "RNA-seq results will have no demographic record to link to: ",
      paste(in_linkage_not_demographics, collapse = ", ")
    )
  }
  if (length(in_demographics_not_linkage) == 0 && length(in_linkage_not_demographics) == 0) {
    cat("Linkage check passed: all participant_id values match between",
        "demographics and linkage files.\n")
  }
}

# ---- WRITE OUTPUT -----------------------------------------------------------

out_path <- file.path(opts$outdir, "cleaned_demographics.tsv")
write_tsv(cleaned, out_path)
cat("Wrote", nrow(cleaned), "rows to", out_path, "\n")

cat("\n--- CONFIRM BEFORE SUBMISSION ---\n")
cat("1. Required fields present:", paste(REQUIRED_FIELDS, collapse = ", "), "\n")
optional_present <- OPTIONAL_FIELDS[OPTIONAL_FIELDS %in% present_target_fields]
if (length(optional_present) > 0) {
  cat("2. Optional fields also included:", paste(optional_present, collapse = ", "), "\n")
} else {
  cat("2. No optional fields (e.g. date_of_birth) found in input -- not included in output.\n")
}
if (!has_linkage) {
  cat("3. No --linkage_file provided -- participant_id/linkage consistency not checked.\n")
}
