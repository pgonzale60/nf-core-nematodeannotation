process HELIXER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://gglyptodon/helixer-docker:latest' :
        'gglyptodon/helixer-docker:latest' }"

    input:
    tuple val(meta), path(fasta)
    val lineage

    output:
    tuple val(meta), path("*_helixer.gff3"), emit: gff3
    tuple val("${task.process}"), val('helixer'), eval("Helixer.py --version 2>&1 | sed -e 's/.*Helixer //; s/ .*//'"), topic: versions, emit: versions_helixer

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def lin      = lineage ?: 'invertebrate'

    """
    export PYTHONNOUSERSITE=1

    Helixer.py \\
        --fasta-path "${fasta}" \\
        --species "${prefix}" \\
        --gff-output-path "${prefix}_helixer.gff3" \\
        --lineage "${lin}" \\
        --temporary-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_helixer.gff3"
    """
}
