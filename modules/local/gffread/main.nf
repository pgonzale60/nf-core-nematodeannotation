process GFFREAD_HELIXER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4':
        'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*.gtf"), emit: gtf
    tuple val("${task.process}"), val('gffread'), val("0.12.7"), emit: versions_gffread, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.helixer"
    """
    gffread -E "${gff3}" -T -o "${prefix}.gtf"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.helixer"
    """
    touch "${prefix}.gtf"
    """
}

process GFFREAD_PROTEIN {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4':
        'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), path(gtf)
    path fasta

    output:
    tuple val(meta), path("*.aa"), emit: aa
    tuple val("${task.process}"), val('gffread'), val("0.12.7"), emit: versions_gffread, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gffread -g "${fasta}" -y "${prefix}.aa" "${gtf}" 2>/dev/null || true
    # Validate that at least one sequence was translated
    if [ ! -s "${prefix}.aa" ]; then
        echo "Error: ${prefix}.aa was not generated or is empty." >&2
        exit 1
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.aa"
    """
}
