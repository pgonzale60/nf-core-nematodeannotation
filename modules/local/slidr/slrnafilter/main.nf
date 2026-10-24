//
// SLIDR_SLRNAFILTER — decide which tail cluster centroids are genuine leader RNAs.
//
// Stage 3 of the tested SLIDR pipeline, split from the BLASTN step because the
// two need different container contents: BLASTN needs the BLAST suite, this
// needs a Python interpreter, and no single container in use here has both.
// Keeping them separate is not just tidiness: the BLASTN step is the expensive
// one (an hour for ~1.4M centroids), and as its own process its result is cached
// independently, so a failure here no longer discards that hour.
//
// THE TEST. A read tail is the 5' overhang the aligner soft-clipped, so for a
// trans-spliced read the tail's 3' end sits exactly at the trans-splice junction,
// just upstream of the leader RNA gene's splice donor. A centroid is therefore a
// leader candidate when its 3' end aligns to the genome (SLIDR_SLRNASITES keeps
// only hits where qend == qlen) AND a splice-donor motif follows within
// --max-overlap bases AND an Sm-binding-site motif follows within --max-rna-len.
//
// Nothing here refers to a known leader sequence. That is the point: this is how
// the workflow answers "what are the leaders in THIS organism", which matters for
// Auanema rhodense, whose canonical SL2 leader is absent from both the assembly
// and the reads.
//
process SLIDR_SLRNAFILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(blast)
    tuple val(meta2), path(centroids)
    tuple val(meta3), path(fasta)
    val donor
    val sm
    val max_rna_len
    val max_overlap

    output:
    tuple val(meta), path("${prefix}.candidate_sl_rnas.fa")  , emit: fasta
    tuple val(meta), path("${prefix}.sl_rna_candidates.tsv") , emit: tsv
    tuple val(meta), path("${prefix}.sl_rna_candidates.gff3"), emit: gff
    tuple val("${task.process}"), val('slidr'), val('0.1.0'), topic: versions, emit: versions_slidr

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    def donor_rg = donor ?: 'GT'
    def sm_rg    = sm ?: '.{20,60}AT{4,6}G'
    def rlength  = max_rna_len ?: 80
    def overlap  = max_overlap ?: 10

    """
    # Tool: slidr_filter_sl_rna.py | Script-SHA: c8eebd9deb16fdec
    slidr_filter_sl_rna.py \\
        --blast ${blast} \\
        --centroids ${centroids} \\
        --fasta ${fasta} \\
        --donor '${donor_rg}' \\
        --sm '${sm_rg}' \\
        --max-rna-len ${rlength} \\
        --max-overlap ${overlap} \\
        --prefix ${prefix} \\
        --outdir . \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.candidate_sl_rnas.fa
    touch ${prefix}.sl_rna_candidates.tsv
    touch ${prefix}.sl_rna_candidates.gff3
    """
}
