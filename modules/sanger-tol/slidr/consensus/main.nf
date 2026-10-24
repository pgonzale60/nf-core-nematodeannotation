//
// SLIDR_CONSENSUS — resolve consensus spliced leaders and trans-splice acceptor sites.
//
// Joins three upstream tables and names leaders DE NOVO in order of use:
//
//     read --(SLIDR_CLUSTERTAILS)--> tail cluster centroid
//     centroid --(SLIDR_SLRNASITES)--> leader RNA gene locus
//
// Nothing about the answer is assumed in advance. The previous version compared
// tails to a hard-coded table of known leaders, which cannot work for an organism
// whose leaders are not in that table -- and Auanema rhodense is one: the
// canonical SL2 leader occurs zero times in its assembly and zero times at read
// 5' ends, while SL1 is abundant.
//
// Retained from the tested implementation: the acceptor-motif test, the rule that
// a trans-spliced gene must sit further than --max-leader-gene-distance bp from
// its own leader gene (otherwise the leader array's own transcripts dominate),
// and the majority-length cluster consensus that separates a leader from
// 5'-truncated fragments of it.
//
process SLIDR_CONSENSUS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sites)
    path readinfo
    path slrnas
    path rrna_gff
    path rfam_gff
    val min_support
    val acceptor
    val max_leader_gene_distance
    val rrna_flank

    output:
    // NOTE: the globs must use '_sl_' (underscore), matching the filenames that
    // slidr_consensus.py actually writes ("${prefix}_sl_acceptor_sites.bed", ...).
    // A dot-separated pattern ("*.sl_acceptor_sites.bed") never matches a
    // "nxAuaRhod1_APS4_transcripts_sl_acceptor_sites.bed" output, which makes the
    // task exit 0 yet fail output collection with MissingFileException.
    tuple val(meta), path("*_sl_acceptor_sites.bed") , emit: bed
    tuple val(meta), path("*_sl_acceptor_sites.gff3"), emit: gff
    tuple val(meta), path("*_sl_summary.tsv")        , emit: summary
    tuple val(meta), path("*_sl_sequences.fa")       , emit: sl_fasta
    tuple val(meta), path("*_sl_site_pool.tsv")      , emit: site_pool
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def min_sup  = min_support ?: 5
    def acc      = acceptor ?: 'AG'
    def ldist    = max_leader_gene_distance ?: 1000
    def r_flank  = rrna_flank ?: 500
    def rrna_opt = rrna_gff ? "--rrna-gff ${rrna_gff} --rrna-flank ${r_flank}" : ""
    def rfam_opt = rfam_gff ? "--rfam-gff ${rfam_gff}" : ""

    """
    # Tool: slidr_consensus.py | Script-SHA: 91ebac3ea90150f9
    slidr_consensus.py \\
        --sites ${sites} \\
        --readinfo ${readinfo} \\
        --slrnas ${slrnas} \\
        ${rrna_opt} \\
        ${rfam_opt} \\
        --acceptor '${acc}' \\
        --min-support ${min_sup} \\
        --max-leader-gene-distance ${ldist} \\
        --prefix "${prefix}" \\
        --outdir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sl_acceptor_sites.bed
    touch ${prefix}_sl_acceptor_sites.gff3
    touch ${prefix}_sl_summary.tsv
    touch ${prefix}_sl_sequences.fa
    touch ${prefix}_sl_site_pool.tsv
    """
}
