process INFERNAL_CMSCANTOGFF3 {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.5--pl5321h7b50bb2_4' :
        'quay.io/biocontainers/infernal:1.1.5--pl5321h7b50bb2_4' }"

    input:
    tuple val(meta), path(tbl)

    output:
    tuple val(meta), path("*.ncrna.gff3")      , emit: gff3
    tuple val(meta), path("*_exclusion.raw.bed"), emit: exclusion_bed
    tuple val("${task.process}"), val('perl'), eval("perl -v | sed -n '2p' | awk '{print \$4}' | tr -d '(),'"), topic: versions, emit: versions_perl

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    perl ${moduleDir}/resources/usr/bin/cmscan_to_gff3.pl \\
        --tbl "${tbl}" \\
        --prefix "${prefix}.ncrna" \\
        --outdir . \\
        ${args}

    mv "${prefix}.ncrna_exclusion.bed" "${prefix}_exclusion.raw.bed"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.ncrna.gff3"
    touch "${prefix}_exclusion.raw.bed"
    """
}
