process INFERNAL_CMSCAN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.5--pl5321h7b50bb2_4' :
        'quay.io/biocontainers/infernal:1.1.5--pl5321h7b50bb2_4' }"

    input:
    tuple val(meta), path(fasta)
    path(cmfile)

    output:
    tuple val(meta), path("*.tbl.gz")     , emit: tbl
    tuple val(meta), path("*.cmscan.log") , emit: log
    tuple val("${task.process}"), val('infernal'), eval("cmscan -h | sed '2!d;s/.*INFERNAL //;s/ .*//'"), topic: versions, emit: versions_infernal

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '--rfam --cut_ga --nohmmonly'
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Tool: cmscan_to_gff3.py | Script-SHA: b3794b7cad738edd
    # Press cmfile if cmpress files not present
    if [ ! -f "${cmfile}.i1m" ]; then
        cmpress "${cmfile}"
    fi

    cmscan \\
        --cpu ${task.cpus} \\
        --tblout "${prefix}.tbl" \\
        ${args} \\
        "${cmfile}" \\
        "${fasta}" > "${prefix}.cmscan.log"

    gzip -c "${prefix}.tbl" > "${prefix}.tbl.gz"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.tbl.gz"
    touch "${prefix}.cmscan.log"
    """
}
