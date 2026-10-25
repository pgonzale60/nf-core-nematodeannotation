process NCRNA_TANDEM_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(ncrna_gff)

    output:
    tuple val(meta), path("*.structural_tandem_arrays.bed"), emit: bed
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Tool: filter_structural_tandem_arrays.py | Script-SHA: 9e13440d175e8de4
    filter_structural_tandem_arrays.py \\
        --gff ${ncrna_gff} \\
        --out ${prefix}.structural_tandem_arrays.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.structural_tandem_arrays.bed
    """
}
