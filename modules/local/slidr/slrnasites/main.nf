//
// SLIDR_SLRNASITES — de-novo discovery of spliced-leader RNA genes.
//
// This is stage 3 of the tested SLIDR pipeline: the step that decides WHAT the
// spliced leaders are in this organism, rather than asserting them from a
// hard-coded table.
//
// WHY THIS MATTERS FOR THIS PROJECT
// ---------------------------------
// The SL1/SL2 sequences supplied in the Agent 2 brief are C. elegans sequences.
// They cannot be assumed for a distantly related nematode, and for Auanema
// rhodense they are demonstrably wrong for SL2: the canonical SL2 leader
// (GGTTTTAACCCAGTTACTCAAG) occurs ZERO times in the curated assembly and zero
// times at read 5' ends across 1M reads, while its 3' 10-mers occur at or below
// random expectation. SL1 by contrast is present (880 copies of the 22-mer,
// 2284 of the 5' 12-mer). A workflow intended to span the genus therefore cannot
// depend on a fixed leader vocabulary; it has to derive one.
//
// HOW IT WORKS (ported verbatim from the tested implementation)
// ------------------------------------------------------------
// A read tail is the 5' overhang the aligner soft-clipped, so for a trans-spliced
// read the tail's 3' end sits exactly at the trans-splice junction, immediately
// upstream of the splice donor of the leader RNA gene. Therefore a tail cluster
// centroid is a genuine leader if and only if:
//
//   1. its 3' end aligns to the genome at 100% identity
//      (enforced upstream by `awk '$8==$14'`, i.e. qend == qlen), and
//   2. a splice-donor motif (GT by default) sits within --max-overlap bases
//      downstream of that alignment, and
//   3. an Sm-binding-site motif (.{20,60}AT{4,6}G by default) follows within
//      --max-rna-len bases.
//
// Low-complexity tails are masked before the search (dustmasker level 10,
// window 10, with -lcase_masking) so that poly-A/poly-T tails cannot be
// mistaken for leaders -- a real risk here, since such tails are abundant.
//
// Outputs feed SLIDR_CONSENSUS, which joins read -> tail cluster -> leader and
// resolves the final acceptor sites.
//
process SLIDR_SLRNASITES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/53/531cf3beb1d695a35e4aae5b80c80ce45d8d46b6128fbe9d07c268cc84956894/data' :
        'community.wave.seqera.io/library/blast_htslib_samtools_seqkit:07b64cc4a2210347' }"

    input:
    tuple val(meta), path(centroids)
    tuple val(meta2), path(fasta)
    val min_reads
    val evalue

    output:
    // This module only MAPS the centroids. Filtering them into candidate leader
    // RNAs runs as SLIDR_SLRNAFILTER in a Python container, because this one
    // (blast + htslib + samtools + seqkit) has no interpreter at all: running the
    // scanner here failed with exit 127 after an hour of BLASTN, and the retry
    // then re-did that whole hour. As a separate process the expensive BLASTN
    // result is cached independently of the cheap filtering step.
    tuple val(meta), path("${prefix}.centroids.blastn.tsv.gz"), emit: blast
    tuple val("${task.process}"), val('blastn'), eval('blastn -version 2>&1 | sed "s/^.*blastn: //; s/ .*$//"'), topic: versions, emit: versions_blast

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def eval_c  = evalue ?: 1
    def minr    = min_reads ?: 5

    """
    mkdir -p blast_db
    makeblastdb -dbtype nucl -in ${fasta} -out blast_db/reference \\
        > blast_db/makeblastdb.log 2>&1

    # Only centroids carrying at least min_reads reads can ever contribute to a
    # reported leader: every read maps to exactly one centroid, a site's support is
    # the number of reads from that centroid, and both the consensus and the
    # reporting rules require >= min_support (the same value) reads. So this filter
    # is lossless for the deliverable rather than a heuristic -- and it is what
    # makes the raised BLASTN E-value affordable: 1,422,818 centroids become
    # ~118,000, and short queries with -evalue 1000 return ~1,700 hits each.
    # The cluster's support is carried in the FASTA header (">centroid_N;size=123").
    awk -v m=$minr '/^>/ { s=\$0; sub(/.*;size=/, "", s); keep = (s+0 >= m); if (keep) print; next } keep' \\
        ${centroids} > centroids_minsupport.fa
    echo "[SLIDR_SLRNASITES] centroids with >= $minr reads: \$(( \$(grep -c '^>' centroids_minsupport.fa) ))"

    # Low-complexity masking, without which homopolymer and other low-complexity
    # centroids (poly-T tails are ~1 % of all tails here) would seed spurious
    # genome-wide hits. dustmasker is the tested path; blastn's own -dust is the
    # fallback so the step degrades rather than silently returning nothing.
    if command -v dustmasker >/dev/null 2>&1; then
        MASK_CMD="dustmasker -in centroids_minsupport.fa -level 10 -window 10 -outfmt fasta"
        MASK_OPT="-lcase_masking"
    else
        echo "[SLIDR_SLRNASITES] dustmasker absent; using blastn -dust 10 10" >&2
        MASK_CMD="cat centroids_minsupport.fa"
        MASK_OPT="-dust 10 10"
    fi

    # Only alignments reaching the 3' end of the centroid (qend == qlen) are
    # informative: the leader's 3' end is what abuts the splice donor, whereas
    # its 5' end may carry sequencing noise.
    \$MASK_CMD \\
    | blastn \\
        -db blast_db/reference \\
        -query - \\
        -word_size 8 \\
        -num_threads ${task.cpus} \\
        -perc_identity 100 \\
        -evalue ${eval_c} \\
        \$MASK_OPT \\
        ${args} \\
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq qlen' \\
        2> ${prefix}.blastn.log \\
    | awk '\$8==\$14' \\
    | gzip -c > ${prefix}.centroids.blastn.tsv.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip -c > ${prefix}.centroids.blastn.tsv.gz
    """
}
