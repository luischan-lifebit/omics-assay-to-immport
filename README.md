# omics-assay-to-immport

Converts raw omics assay output into ImmPort's submission templates.

Part of the assay to ImmPort to OMOP pipeline ([Miro board](https://miro.com/app/board/uXjVHxn9PL8=/)). This repo covers the "make into ImmPort format" step; output feeds into [`immport_to_omop`](https://github.com/lifebit-ai/immport_to_omop).

## `bin/salmon_to_immport_rnaseq.R`

Converts merged Salmon RNA-seq TPM matrices into ImmPort's `RNA_SEQ_Results` template (schema 3.37).

**Input:** `salmon.merged.gene_tpm.tsv`, `salmon.merged.transcript_tpm.tsv`
**Output:** `RNA_SEQ_Results_gene.tsv`, `RNA_SEQ_Results_transcript.tsv`

One row per sample x gene/transcript, in ImmPort's required format. Controlled vocab (Repository Name, Transcript Type, Result Unit) is pulled directly from ImmPort's official template, not guessed.

## Running as a Nextflow pipeline

Requires Nextflow >= 22.10 and either Conda or Docker.

```bash
# with conda
nextflow run main.nf \
    --gene_tpm salmon.merged.gene_tpm.tsv \
    --transcript_tpm salmon.merged.transcript_tpm.tsv \
    --outdir results \
    -profile conda

# with docker
nextflow run main.nf \
    --gene_tpm salmon.merged.gene_tpm.tsv \
    --transcript_tpm salmon.merged.transcript_tpm.tsv \
    --outdir results \
    -profile docker
```

### Parameters

| Param | Default | Description |
|---|---|---|
| `--gene_tpm` | *(required)* | Path to `salmon.merged.gene_tpm.tsv` |
| `--transcript_tpm` | *(required)* | Path to `salmon.merged.transcript_tpm.tsv` |
| `--outdir` | `results` | Output directory |
| `--repository_name` | `Ensembl` | Repository Name value (still a guess, see open items above) |
| `--transcript_type` | `mRNA` | Transcript Type Reported, applied to all rows |
| `--result_unit` | `TPM` | Result Unit Reported |

### Output

- `results/RNA_SEQ_Results_gene.tsv`
- `results/RNA_SEQ_Results_transcript.tsv`
- `results/pipeline_info/` — execution report and trace, for reproducibility records

## Running the R script directly (no Nextflow)

```bash
Rscript bin/salmon_to_immport_rnaseq.R \
    salmon.merged.gene_tpm.tsv \
    salmon.merged.transcript_tpm.tsv \
    --outdir results
```

Requires: `dplyr`, `tidyr`, `readr`, `tibble`, `optparse`

```r
install.packages(c("dplyr", "tidyr", "readr", "tibble", "optparse"))
```

Run output through the [ImmPort Validator](https://docs.immport.org/datasubmission/datapackagevalidator/) before any real upload.

## Repo structure

```
omics-assay-to-immport/
├── README.md
├── main.nf              # pipeline definition
├── nextflow.config       # params, conda/docker profiles
├── environment.yml       # R + package versions for the conda profile
└── bin/
    └── salmon_to_immport_rnaseq.R   # conversion script (Nextflow auto-adds bin/ to PATH)
```
