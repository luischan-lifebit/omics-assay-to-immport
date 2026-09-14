# omics-assay-to-immport

Converts raw omics assay output into ImmPort's submission templates, and cleans demographic files for OMOP compatibility.

Part of the assay to ImmPort to OMOP pipeline ([Miro board](https://miro.com/app/board/uXjVHxn9PL8=/)). This repo covers the "make into ImmPort format" step; output feeds into [`immport_to_omop`](https://github.com/lifebit-ai/immport_to_omop).

This pipeline has two independent capabilities, either or both runnable in a single invocation:

1. **RNA-seq conversion** (`--gene_tpm` / `--transcript_tpm`) — converts Salmon TPM matrices into ImmPort's `RNA_SEQ_Results` format.
2. **Demographics cleaning** (`--demographics_file`) — reformats a demographic file's column headers to match what `immport_to_omop`'s field map expects.

## `bin/salmon_to_immport_rnaseq.R`

Converts merged Salmon RNA-seq TPM matrices into ImmPort's `RNA_SEQ_Results` template (schema 3.37).

**Input:** `salmon.merged.gene_tpm.tsv`, `salmon.merged.transcript_tpm.tsv`
**Output:** `RNA_SEQ_Results_gene.tsv`, `RNA_SEQ_Results_transcript.tsv`

One row per sample x gene/transcript. Column headers are lowercase_with_underscores (e.g. `expsample_id`), matching `source_column_name` in Lifebit's maintained OMOP field map templates directly — no manual renaming needed before OMOP ETL. Controlled vocab (Repository Name, Transcript Type, Result Unit) is pulled directly from ImmPort's official template, not guessed.

Optional participant linkage: when `--linkage_file` is provided (a `source_person_id,sample_id` CSV, per Lifebit's `lifebit_omics_linkage` contract), a `participant_id` column is added, joined from `sample_id` to `source_person_id`. Samples with no matching linkage entry are flagged with a warning rather than silently dropped.

## `bin/clean_demographics.R`

Cleans a user-provided demographic file's column headers so it's directly compatible with Data Factory OMOP ETL. Does not generate or infer missing demographic data — only reformats what's provided.

**Input:** any CSV/TSV demographic file with the required fields, in any reasonably-labelled raw form (e.g. `Participant ID`, `participant id`, `PARTICIPANT_ID` all match).
**Output:** `cleaned_demographics.tsv`

Required fields: `participant_id`, `gender`, `year_of_birth`, `race`, `ethnicity` — confirmed directly against the OMOP person table's actual Required validation. `date_of_birth` is included if present but is not required (`year_of_birth` is the required date field, not `date_of_birth`/`birth_datetime`).

If a required field is missing, the process fails with a clear message naming exactly which column(s) are missing, rather than proceeding with incomplete data.

Column matching is exact-label only (case/whitespace/underscore-insensitive) — the same convention as the RNA-seq script. Synonyms (e.g. `Sex` for `gender`, `DOB` for `year_of_birth`) are not auto-matched.

When `--linkage_file` is also provided, this script checks whether the demographic file's `participant_id` values overlap with the linkage file's `source_person_id` values, and **warns** (does not fail) on any mismatch in either direction.

## `bin/rename_headers.R`

Generic column-header renamer, used to produce ImmPort's canonical submission format (Title Case headers) from the RNA-seq script's lowercase_with_underscores output. Takes any TSV and a mapping CSV (`old_name,new_name`) and writes a renamed copy — it knows nothing about RNA-seq, ImmPort, or OMOP specifically, so the same script is reusable for other assay types by supplying a different mapping file, without touching code.

Any column in the input not listed in the mapping file is passed through unchanged, rather than dropped.

## Running as a Nextflow pipeline

Requires Nextflow >= 22.10 and either Conda or Docker.

```bash
# RNA-seq conversion only
nextflow run main.nf \
    --gene_tpm salmon.merged.gene_tpm.tsv \
    --transcript_tpm salmon.merged.transcript_tpm.tsv \
    --outdir results \
    -profile docker

# Demographics cleaning only
nextflow run main.nf \
    --demographics_file demographics.csv \
    --outdir results \
    -profile docker

# Both together, with participant linkage
nextflow run main.nf \
    --gene_tpm salmon.merged.gene_tpm.tsv \
    --transcript_tpm salmon.merged.transcript_tpm.tsv \
    --demographics_file demographics.csv \
    --linkage_file lifebit_omics_linkage.csv \
    --outdir results \
    -profile docker
```

### Parameters

| Param | Default | Description |
|---|---|---|
| `--gene_tpm` | *(none)* | Path to `salmon.merged.gene_tpm.tsv`. Required together with `--transcript_tpm` to run RNA-seq conversion. |
| `--transcript_tpm` | *(none)* | Path to `salmon.merged.transcript_tpm.tsv`. Required together with `--gene_tpm`. |
| `--demographics_file` | *(none)* | Path to a demographic file. Runs demographics cleaning independently of the RNA-seq inputs. |
| `--outdir` | `results` | Output directory |
| `--repository_name` | `Ensembl` | Repository Name value for RNA-seq output (still a guess, see open items below) |
| `--transcript_type` | `mRNA` | Transcript Type Reported, applied to all RNA-seq rows |
| `--result_unit` | `TPM` | Result Unit Reported for RNA-seq output |
| `--linkage_file` | *(none)* | Optional `source_person_id,sample_id` CSV, used by both RNA-seq conversion (adds `participant_id`) and demographics cleaning (validates `participant_id` consistency) |
| `--header_map` | `assets/immport_header_map.csv` | Mapping file used to produce the ImmPort submission format from RNA-seq output |

At least one of `--demographics_file` or (`--gene_tpm` + `--transcript_tpm`) is required. `--gene_tpm` and `--transcript_tpm` must be provided together.

### Output
```
results/
├── RNA_SEQ_Results_gene.tsv # CloudOS/OMOP-ready, if RNA-seq inputs given
├── RNA_SEQ_Results_transcript.tsv # CloudOS/OMOP-ready
├── cleaned_demographics.tsv # CloudOS/OMOP-ready, if --demographics_file given
├── immport_original_format/ # RNA-seq only -- same data, ImmPort's canonical Title Case headers
│ ├── RNA_SEQ_Results_gene.tsv
│ └── RNA_SEQ_Results_transcript.tsv
└── pipeline_info/ # execution report and trace, for reproducibility records
```



Run the `immport_original_format/` output through the [ImmPort Validator](https://docs.immport.org/datasubmission/datapackagevalidator/) before any real ImmPort submission.

## Running the R scripts directly (no Nextflow)

```bash
Rscript bin/salmon_to_immport_rnaseq.R \
    salmon.merged.gene_tpm.tsv \
    salmon.merged.transcript_tpm.tsv \
    --outdir results

Rscript bin/clean_demographics.R \
    demographics.csv \
    --outdir results

Rscript bin/rename_headers.R \
    results/RNA_SEQ_Results_gene.tsv \
    assets/immport_header_map.csv \
    results/immport_original_format/RNA_SEQ_Results_gene.tsv
```

Requires: `dplyr`, `tidyr`, `readr`, `tibble`, `stringr`

```r
install.packages(c("dplyr", "tidyr", "readr", "tibble", "stringr"))
```

## Testing

Two demographic test fixtures are provided under `test_data/`:

- `demographics_dirty_example.csv` — badly-formatted but complete headers; should clean successfully.
- `demographics_missing_column.csv` — missing required fields; should fail with a clear message naming what's missing.

```bash
Rscript bin/clean_demographics.R test_data/demographics_dirty_example.csv --outdir test_output
Rscript bin/clean_demographics.R test_data/demographics_missing_column.csv --outdir test_output
```

These are fixtures for validating the cleaning logic only — not synthetic data generated by the pipeline itself, and not used anywhere in the real pipeline.

## Repo structure

```
omics-assay-to-immport/
├── README.md
├── main.nf # pipeline definition
├── nextflow.config # params, conda/docker profiles
├── environment.yml # R + package versions for the conda profile
├── assets/
│ ├── immport_header_map.csv # old_name,new_name mapping used by rename_headers.R
│ └── NO_LINKAGE_FILE # placeholder when --linkage_file isn't provided
├── test_data/
│ ├── demographics_dirty_example.csv
│ └── demographics_missing_column.csv
└── bin/
├── salmon_to_immport_rnaseq.R # RNA-seq TPM matrices -> ImmPort format
├── clean_demographics.R # demographic file header cleaning
└── rename_headers.R # generic header renamer (Nextflow auto-adds bin/ to PATH)
```
