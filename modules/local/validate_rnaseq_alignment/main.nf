process VALIDATE_RNASEQ_ALIGNMENT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(alignment), path(index)
    tuple val(meta_fai), path(fai)

    output:
    tuple val(meta), path(alignment), path(index), env('IS_HISAT2'), env('IS_COORDINATE_VALID'), emit: evaluation
    tuple val(meta), path(alignment), path(index), env('IS_SPLICED'), env('IS_VALID'), emit: spliced_evaluation, optional: true
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def reference_opt = ""
    def fasta_candidate = file("${fai.toString().replaceAll(/\\.fai\$/, '')}")
    if (fasta_candidate.exists()) {
        reference_opt = "--reference ${fasta_candidate}"
    }
    """
    # Tool: validate_rnaseq_alignment.py | Script-SHA: 4e43f1e1befbe9e6
    validate_rnaseq_alignment.py \\
        --alignment "${alignment}" \\
        --fai "${fai}" \\
        ${reference_opt} \\
        --env-out validation.env

    source validation.env
    """

    stub:
    """
    IS_HISAT2="true"
    IS_COORDINATE_VALID="true"
    IS_SPLICED="true"
    IS_VALID="true"
    """
}
