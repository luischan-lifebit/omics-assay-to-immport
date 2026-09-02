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
 *
 * Optional participant linkage (adds a Participant ID column, joined from
 * sample_id -> source_person_id, per Lifebit's lifebit_omics_linkage
 * contract):
 *   nextflow run main.nf \
 *       --gene_tpm salmon.merged.gene_tpm.tsv \
 *       --transcript_tpm salmon.merged.transcript_tpm.tsv \
 *       --linkage_file lifebit_omics_linkage.csv \
 *       --outdir results
 */

params.gene_tpm         = null
params.transcript_tpm   = null
params.outdir           = "results"
params.repository_name  = "Ensembl"
params.transcript_type  = "mRNA"
params.result_unit      = "TPM"
params.linkage_file     = null   // optional: source_person_id,sample_id CSV

if (!params.gene_tpm || !params.transcript_tpm) {
    error """
    Missing required input.
    Usage: nextflow run main.nf --gene_tpm <path> --transcript_tpm <path> [--outdir results] [--linkage_file <path>]
    """
}

process SALMON_TO_IMMPORT_RNASEQ {
    tag "rnaseq_to_immport"
    publishDir params.outdir, mode: 'copy'

    input:
    path gene_tpm
    path transcript_tpm
    path linkage_file

    output:
    path "RNA_SEQ_Results_gene.tsv",       emit: gene_results
    path "RNA_SEQ_Results_transcript.tsv", emit: transcript_results

    script:
    def linkage_arg = params.linkage_file ? "--linkage_file ${linkage_file}" : ""
    """
    salmon_to_immport_rnaseq.R \\
        ${gene_tpm} \\
        ${transcript_tpm} \\
        --repository_name "${params.repository_name}" \\
        --transcript_type "${params.transcript_type}" \\
        --result_unit "${params.result_unit}" \\
        --outdir . \\
        ${linkage_arg}
    """
}

workflow {
    gene_tpm_ch       = Channel.fromPath(params.gene_tpm, checkIfExists: true)
    transcript_tpm_ch = Channel.fromPath(params.transcript_tpm, checkIfExists: true)

    // Use a real placeholder file when no linkage file is given, since the
    // process declares a fixed 'path linkage_file' input (avoids the
    // "optional:" syntax that broke on this Nextflow version).
    if (params.linkage_file) {
        linkage_ch = Channel.fromPath(params.linkage_file, checkIfExists: true)
    } else {
        linkage_ch = Channel.fromPath("${projectDir}/assets/NO_LINKAGE_FILE")
    }

    SALMON_TO_IMMPORT_RNASEQ(gene_tpm_ch, transcript_tpm_ch, linkage_ch)
}

workflow.onComplete {
    log.info """
    Pipeline complete: ${workflow.success ? 'OK' : 'FAILED'}
    Output directory : ${params.outdir}
    Linkage file used: ${params.linkage_file ?: '(none -- no Participant ID column added)'}

    Reminders (still open, per README):
      1. Repository Name  = '${params.repository_name}' -- confirm correct for this data's ID system.
      2. Transcript Type   = '${params.transcript_type}' applied to ALL rows -- needs real biotype data.
      ${params.linkage_file ? "3. Confirm linkage file's source_person_id values match person.person_source_value in the target OMOP schema." : "3. No linkage file provided -- add one if this feeds into OMOP ingestion requiring participant linkage."}
    Run output through the ImmPort Validator before any real upload.
    """
}
