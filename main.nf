#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * omics-assay-to-immport
 *
 * Two independent capabilities, either or both runnable in a single
 * invocation:
 *
 *   1. RNA-seq conversion (--gene_tpm / --transcript_tpm)
 *        SALMON_TO_IMMPORT_RNASEQ  -- TPM matrices -> ImmPort long format
 *        IMMPORT_DATA_MODEL_TO_OMOP_COMPATIBLE_FORMAT  -- generic header rename -> ImmPort
 *                                     submission format (Title Case)
 *
 *   2. Demographics cleaning (--demographics_file)
 *        CLEAN_DEMOGRAPHICS -- reformats a user-provided demographic file's
 *        headers to match what immport_to_omop's field map expects. Does
 *        not generate or infer missing data; fails with a clear message if
 *        a required column is missing.
 *
 * --linkage_file is shared across both: it adds a participant_id column to
 * the RNA-seq output, and is used by CLEAN_DEMOGRAPHICS to check that
 * demographic participant_id values match the linkage file's
 * source_person_id values.
 *
 * Usage:
 *   nextflow run main.nf --gene_tpm <path> --transcript_tpm <path> --outdir results
 *   nextflow run main.nf --demographics_file <path> --outdir results
 *   nextflow run main.nf --gene_tpm <path> --transcript_tpm <path> \
 *       --demographics_file <path> --linkage_file <path> --outdir results
 */

params.gene_tpm          = null
params.transcript_tpm    = null
params.demographics_file = null
params.outdir            = "results"
params.repository_name   = "Ensembl"
params.transcript_type   = "mRNA"
params.result_unit       = "TPM"
params.linkage_file      = null   // optional: source_person_id,sample_id CSV
params.header_map        = "${projectDir}/assets/immport_header_map.csv"

def run_rnaseq       = params.gene_tpm && params.transcript_tpm
def run_demographics  = params.demographics_file as boolean

if (!run_rnaseq && !run_demographics) {
    error """
    Missing required input.
    Provide either --gene_tpm + --transcript_tpm, or --demographics_file (or both).
    """
}
if ((params.gene_tpm && !params.transcript_tpm) || (!params.gene_tpm && params.transcript_tpm)) {
    error "Both --gene_tpm and --transcript_tpm are required together."
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

process RENAME_TO_IMMPORT_FORMAT {
    tag "${input_tsv.baseName}"
    publishDir "${params.outdir}/immport_original_format", mode: 'copy'

    input:
    path input_tsv
    path header_map

    output:
    path "${input_tsv.name}", emit: renamed

    script:
    """
    rename_headers.R ${input_tsv} ${header_map} ${input_tsv.name}
    """
}

process CLEAN_DEMOGRAPHICS {
    tag "clean_demographics"
    publishDir params.outdir, mode: 'copy'

    input:
    path demographics_file
    path linkage_file

    output:
    path "cleaned_demographics.tsv", emit: cleaned_demographics

    script:
    def linkage_arg = params.linkage_file ? "--linkage_file ${linkage_file}" : ""
    """
    clean_demographics.R \\
        ${demographics_file} \\
        --outdir . \\
        ${linkage_arg}
    """
}

workflow {
    // shared linkage channel -- real placeholder when not provided, since
    // both processes declare a fixed 'path linkage_file' input
    if (params.linkage_file) {
        linkage_ch = Channel.fromPath(params.linkage_file, checkIfExists: true)
    } else {
        linkage_ch = Channel.fromPath("${projectDir}/assets/NO_LINKAGE_FILE")
    }

    if (run_rnaseq) {
        gene_tpm_ch       = Channel.fromPath(params.gene_tpm, checkIfExists: true)
        transcript_tpm_ch = Channel.fromPath(params.transcript_tpm, checkIfExists: true)
        header_map_ch     = Channel.fromPath(params.header_map, checkIfExists: true)

        SALMON_TO_IMMPORT_RNASEQ(gene_tpm_ch, transcript_tpm_ch, linkage_ch)

        rnaseq_outputs_ch = SALMON_TO_IMMPORT_RNASEQ.out.gene_results
            .mix(SALMON_TO_IMMPORT_RNASEQ.out.transcript_results)

        RENAME_TO_IMMPORT_FORMAT(rnaseq_outputs_ch, header_map_ch)
    }

    if (run_demographics) {
        demographics_ch = Channel.fromPath(params.demographics_file, checkIfExists: true)
        CLEAN_DEMOGRAPHICS(demographics_ch, linkage_ch)
    }
}

workflow.onComplete {
    log.info """
    Pipeline complete: ${workflow.success ? 'OK' : 'FAILED'}
    Output directory : ${params.outdir}
    RNA-seq conversion run   : ${run_rnaseq}
    Demographics cleaning run: ${run_demographics}
    Linkage file used: ${params.linkage_file ?: '(none)'}

    ${run_rnaseq ? """RNA-seq reminders:
      - Repository Name = '${params.repository_name}' -- confirm correct for this data's ID system.
      - Transcript Type  = '${params.transcript_type}' applied to ALL rows -- needs real biotype data.
      Run the immport_original_format/ output through the ImmPort Validator before any real upload.""" : ""}
    ${run_demographics ? """Demographics reminders:
      - Required fields: participant_id, gender, year_of_birth, race, ethnicity.
      - date_of_birth -> birth_datetime is included if present, but not required.
      ${params.linkage_file ? "- Check the log above for any participant_id/linkage mismatch warnings." : "- No linkage file provided -- participant_id/linkage consistency was not checked."}""" : ""}
    """
}
