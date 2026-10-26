process SATELLITEFILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(dat_files)
    val min_period
    val min_copies
    val min_identity

    output:
    tuple val(meta), path("*_filtered_units.fa"), emit: fasta
    tuple val(meta), path("*_trf_summary.tsv")  , emit: summary, optional: true
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def period   = min_period ?: 50
    def copies   = min_copies ?: 100
    def identity = min_identity ?: 85

    """
    # Tool: filter_trf_satellites.py | Script-SHA: 8cac59fe405af914
    python3 ${moduleDir}/resources/usr/bin/filter_trf_satellites.py \\
        --dat ${dat_files} \\
        --min-period ${period} \\
        --min-copies ${copies} \\
        --min-identity ${identity} \\
        --prefix "${prefix}" \\
        --outdir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_filtered_units.fa
    touch ${prefix}_trf_summary.tsv
    """
}
