#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * omics-assay-to-immport
 * Converts Salmon RNA-seq TPM output into ImmPort's RNA_SEQ_Results format.
 *
 * Two processes:
 *   1. SALMON_TO_IMMPORT_RNASEQ -- the science: TPM matrices -> ImmPort long
 *      format, lowercase_with_underscores headers (CloudOS/OMOP-ready).
 *   2. RENAME_TO_IMMPORT_FORMAT -- generic header renaming: takes process 1's
 *      output + a header-map asset, produces the canonical ImmPort
 *      submission format (Title Case headers) alongside it. This step knows
 *      nothing about RNA-seq specifically, so it's reusable as-is for other
 *      assay types (proteomics, demographics) by swapping the header-map
 *      asset -- splitting format concerns from conversion logic, at the
 *      workflow level rather than inside each conversion script.
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
params.header_map       = "${projectDir}/assets/immport_header_map.csv"

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

workflow {
    gene_tpm_ch       = Channel.fromPath(params.gene_tpm, checkIfExists: true)
    transcript_tpm_ch = Channel.fromPath(params.transcript_tpm, checkIfExists: true)
    header_map_ch     = Channel.fromPath(params.header_map, checkIfExists: true)

    // Use a real placeholder file when no linkage file is given, since the
    // process declares a fixed 'path linkage_file' input (avoids the
    // "optional:" syntax that broke on this Nextflow version).
    if (params.linkage_file) {
        linkage_ch = Channel.fromPath(params.linkage_file, checkIfExists: true)
    } else {
        linkage_ch = Channel.fromPath("${projectDir}/assets/NO_LINKAGE_FILE")
    }

    SALMON_TO_IMMPORT_RNASEQ(gene_tpm_ch, transcript_tpm_ch, linkage_ch)

    // Run the generic renamer on both outputs, each paired with the same
    // header-map asset.
    rnaseq_outputs_ch = SALMON_TO_IMMPORT_RNASEQ.out.gene_results
        .mix(SALMON_TO_IMMPORT_RNASEQ.out.transcript_results)

    RENAME_TO_IMMPORT_FORMAT(rnaseq_outputs_ch, header_map_ch)
}

workflow.onComplete {
    log.info """
    Pipeline complete: ${workflow.success ? 'OK' : 'FAILED'}
    Output directory : ${params.outdir}
      - CloudOS/OMOP format:            ${params.outdir}/RNA_SEQ_Results_{gene,transcript}.tsv
      - ImmPort submission format:      ${params.outdir}/immport_original_format/RNA_SEQ_Results_{gene,transcript}.tsv
    Linkage file used: ${params.linkage_file ?: '(none -- no Participant ID column added)'}

    Reminders (still open, per README):
      1. Repository Name  = '${params.repository_name}' -- confirm correct for this data's ID system.
      2. Transcript Type   = '${params.transcript_type}' applied to ALL rows -- needs real biotype data.
      ${params.linkage_file ? "3. Confirm linkage file's source_person_id values match person.person_source_value in the target OMOP schema." : "3. No linkage file provided -- add one if this feeds into OMOP ingestion requiring participant linkage."}
    Run the immport_original_format/ output through the ImmPort Validator before any real upload.
    """
}
