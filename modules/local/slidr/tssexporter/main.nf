//
// SLIDR_TSSEXPORTER — export empirical trans-splicing TSS hints from the SLIDR
// locus-level pool (*_sl_site_pool.tsv written by SLIDR_CONSENSUS).
//
// Thresholds apply to reads pooled per genomic site, not per tail-cluster
// centroid. Output is GFF3, suitable as hint input for gene callers.
//
process SLIDR_TSSEXPORTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(pool)
    val min_reads
    val min_libraries

    output:
    tuple val(meta), path("*.tss.gff3")       , emit: tss_gff
    tuple val(meta), path("*.tss_metrics.tsv"), emit: metrics
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def min_r    = min_reads ?: 3
    def min_l    = min_libraries ?: 1

    """
    # Tool: slidr_tss_exporter.py | Script-SHA: fe68134bed185d4c
    slidr_tss_exporter.py \\
        --pool ${pool} \\
        --prefix "${prefix}" \\
        --min-reads ${min_r} \\
        --min-libraries ${min_l} \\
        --outdir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tss.gff3
    touch ${prefix}.tss_metrics.tsv
    """
}
