process SLOPPR_OPERONS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bed)
    path gff
    val max_intercistronic_dist

    output:
    // Globs must use '_predicted_operons' / '_operons_summary' (underscore) to
    // match what sloppr_operons.py writes in the same directory:
    //   "${prefix}_predicted_operons.gff3"
    //   "${prefix}_operons_summary.tsv"
    // Dot-separated patterns ("*.predicted_operon*.gff3") never match those
    // names, which makes the task exit 0, leave correct files behind, and then
    // fail output collection with MissingFileException. This is the third
    // occurrence of that defect in this workflow's modules; the mismatch was
    // found by a -stub-run of the SPLICED_LEADER_OPERONS subworkflow.
    tuple val(meta), path("*_predicted_operons.gff3"), emit: operons_gff
    tuple val(meta), path("*_operons_summary.tsv")   , emit: summary, optional: true
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def gff_arg  = gff ? "--gff ${gff}" : ""
    def max_dist = max_intercistronic_dist ?: 1000

    """
    # Tool: sloppr_operons.py | Script-SHA: b8737da9eb2cd1ae
    sloppr_operons.py \\
        --bed "${bed}" \\
        ${gff_arg} \\
        --max-intercistronic-dist ${max_dist} \\
        --prefix "${prefix}" \\
        --outdir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predicted_operons.gff3
    touch ${prefix}_operons_summary.tsv
    """
}
