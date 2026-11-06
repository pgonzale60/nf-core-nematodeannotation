process SLIDR_EXTRACTTAILS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    val min_len

    output:
    tuple val(meta), path("*.tails.fa.gz")          , emit: tails
    tuple val(meta), path("*.acceptor_sites.tsv.gz"), emit: sites
    tuple val(meta), path("*.tail_stats.json")      , emit: stats
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stranded = task.ext.strandedness ?: "0"

    """
    # Tool: slidr_extract_tails.py | Script-SHA: b2f5da18cba26d5d
    slidr_extract_tails.py \\
        --bam "${bam}" \\
        --fasta "${fasta}" \\
        --prefix "${prefix}" \\
        --min-tail-len ${min_len} \\
        --strandedness "${stranded}" \\
        --outdir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip -c > ${prefix}.tails.fa.gz
    printf "read_id\\tchrom\\tjunction_pos_0based\\tjunction_pos_1based\\tstrand\\ttail_type\\tclip_len\\ttail_seq\\tgenomic_dinuc\\tis_canonical_ag\\n" | gzip -c > ${prefix}.acceptor_sites.tsv.gz
    echo '{"total_processed_alignments":0,"extracted_sl_tails":0,"canonical_ag_tails":0,"canonical_ag_percentage":0.0,"unique_junction_sites":0,"min_tail_length":8,"strandedness":"0"}' > ${prefix}.tail_stats.json
    """
}
