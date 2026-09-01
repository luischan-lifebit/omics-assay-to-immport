#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * omics-assay-to-immport
 * Converts Salmon RNA-seq TPM output into ImmPort's RNA_SEQ_Results format.
 *
 * Usage:
 *   nextflow run main.nf \
 *       --gene_tpm salmon.merged.gene_tpm.tsv \
 *       --transcript_tpm salmon.merged.transcript_tpm.tsv \
 *       --outdir results
 */

params.gene_tpm         = null
params.transcript_tpm   = null
params.outdir           = "results"
params.repository_name  = "Ensembl"
params.transcript_type  = "mRNA"
params.result_unit      = "TPM"

if (!params.gene_tpm || !params.transcript_tpm) {
    error """
    Missing required input.
    Usage: nextflow run main.nf --gene_tpm <path> --transcript_tpm <path> [--outdir results]
    """
}

process SALMON_TO_IMMPORT_RNASEQ {
    tag "rnaseq_to_immport"
    publishDir params.outdir, mode: 'copy'

    input:
    path gene_tpm
    path transcript_tpm

    output:
    path "RNA_SEQ_Results_gene.tsv",       emit: gene_results
    path "RNA_SEQ_Results_transcript.tsv", emit: transcript_results

    script:
    """
    salmon_to_immport_rnaseq.R \\
        ${gene_tpm} \\
        ${transcript_tpm} \\
        --repository_name "${params.repository_name}" \\
        --transcript_type "${params.transcript_type}" \\
        --result_unit "${params.result_unit}" \\
        --outdir .
    """
}

workflow {
    gene_tpm_ch       = Channel.fromPath(params.gene_tpm, checkIfExists: true)
    transcript_tpm_ch = Channel.fromPath(params.transcript_tpm, checkIfExists: true)

    SALMON_TO_IMMPORT_RNASEQ(gene_tpm_ch, transcript_tpm_ch)
}

workflow.onComplete {
    log.info """
    Pipeline complete: ${workflow.success ? 'OK' : 'FAILED'}
    Output directory : ${params.outdir}

    Reminders (still open, per README):
      1. Repository Name  = '${params.repository_name}' -- confirm correct for this data's ID system.
      2. Transcript Type   = '${params.transcript_type}' applied to ALL rows -- needs real biotype data.
      3. Expsample ID      = raw sample names -- swap for real ImmPort IDs once available.
    Run output through the ImmPort Validator before any real upload.
    """
}
