process SATELLITEPROFILER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bed), path(centroids)
    val max_gap

    output:
    tuple val(meta), path("*_satellite_arrays.gff3")  , emit: arrays_gff3
    tuple val(meta), path("*_satellite_monomers.gff3"), emit: monomers_gff3
    tuple val(meta), path("*_satellite_families.tsv") , emit: families_tsv
    tuple val(meta), path("*_satellite_monomers.fasta"), emit: library_fasta
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def gap_val = max_gap ?: 100

    """
    # Tool: satellite_profiler.py | Script-SHA: 3d80096bbae9f61f
    python3 ${moduleDir}/resources/usr/bin/satellite_profiler.py \\
        --bed "${bed}" \\
        --centroids "${centroids}" \\
        --max-gap ${gap_val} \\
        --prefix "${prefix}" \\
        --outdir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_satellite_arrays.gff3
    touch ${prefix}_satellite_monomers.gff3
    touch ${prefix}_satellite_families.tsv
    touch ${prefix}_satellite_monomers.fasta
    """
}
