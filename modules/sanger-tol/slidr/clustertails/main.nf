//
// SLIDR_CLUSTERTAILS — dereplicate and 3'-anchor-cluster read tails.
//
// Stage 2 of the tested SLIDR pipeline, replacing
//
//     vsearch --derep_fulllength | vsearch --cluster_fast --id 1.0 --iddef 0 --rightjust
//
// with a dependency-free implementation (see slidr_cluster_tails.py for the full
// rationale). Two properties are essential and easy to get wrong:
//
//   * clustering must be 3'-anchored. A tail's 3' end sits at the trans-splice
//     junction, so tails belong together when their 3' ends agree, no matter how
//     much of the leader's 5' end survived. In vsearch terms that is
//     --rightjust combined with --iddef 0 (identity over the shorter sequence),
//     which is exactly "the shorter tail is an exact suffix of the longer".
//
//   * tails must be pooled ACROSS libraries before clustering, because a site
//     supported by 3 reads in each of 16 libraries is a real site but no single
//     library would cluster it.
//
// Input is the collected set of per-library *.tails.fa.gz files. Output is the
// cluster centroids (for SL RNA discovery), a per-cluster summary, and the
// read -> centroid map that SLIDR_CONSENSUS needs to attribute every acceptor
// site to a discovered leader.
//
process SLIDR_CLUSTERTAILS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(tails)
    val min_len

    output:
    tuple val(meta), path("${prefix}_centroids.fa")       , emit: centroids
    tuple val(meta), path("${prefix}_cluster_summary.tsv"), emit: summary
    tuple val(meta), path("${prefix}_readinfo.tsv.gz")    , emit: readinfo
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def minl = min_len ?: 8

    """
    # Tool: slidr_cluster_tails.py | Script-SHA: 5f21be064fa1d26e
    slidr_cluster_tails.py \\
        --tails ${tails} \\
        --min-len ${minl} \\
        --prefix ${prefix} \\
        --outdir . \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_centroids.fa
    printf 'centroid_id\tsequence\treads\tunique_tails\tlength\n' > ${prefix}_cluster_summary.tsv
    echo "" | gzip -c > ${prefix}_readinfo.tsv.gz
    """
}
