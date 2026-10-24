process ARCHIVE_SUBOPTIMAL_MODELS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:22.04'
        : 'quay.io/biocontainers/ubuntu:22.04'}"

    input:
    tuple val(meta), path(gtfs)

    output:
    tuple val(meta), path("${prefix}.suboptimal_gene_models.tar.gz"), emit: archive
    tuple val("${task.process}"), val('tar'), eval("tar --version | head -n 1 | awk '{print \$4}'"), topic: versions, emit: versions_tar

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p suboptimal_gene_models
    for gtf in ${gtfs}; do
        cp -L "\${gtf}" suboptimal_gene_models/
    done

    tar -czf "${prefix}.suboptimal_gene_models.tar.gz" suboptimal_gene_models/
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.suboptimal_gene_models.tar.gz"
    """
}
