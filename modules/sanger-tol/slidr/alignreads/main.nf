//
// HISAT2_ALIGN_SOFTCLIP — soft-clip-permissive spliced read alignment.
// Serves as the primary RNA-seq aligner feeding spliced leader mining and BRAKER3.
// Formerly named SLIDR_ALIGNREADS.
//
// WHY THIS EXISTS RATHER THAN HISAT2_ALIGN
// ----------------------------------------
// SLIDR/SLOPPR recover trans-splice acceptor sites from the 5' sequences that the
// aligner soft-clips off RNA-seq reads. That only works if the aligner is told to
// prefer clipping over mismatching, which needs the scoring options SLIDR
// prescribes (Wenzel et al. 2021, README "soft-clip alignment"):
//
//     hisat2 --sp 1,0 --mp 3,1 -k 5 --score-min L,5,-0.4*<scale>
//
// --sp 1,0 makes the minimum score equal to the maximum number of soft-clipped
// bases, so --score-min L,5,-0.4 caps clipping at ~35 bp for a 100 bp read,
// which comfortably fits a nematode SL.
//
// The stock modules/nf-core/hisat2/align module is NOT usable here for two
// independent reasons:
//   1. it adds --no-mixed --no-discordant for paired input, discarding any pair
//      whose mate failed to map, and
//   2. it pipes through `samtools view -bS -F 4 -F 8 -F 256`, which drops
//      secondary alignments (0x100).
// The tested SLIDR pipeline does neither: slidr_library.sh keeps secondary
// alignments on purpose ("we want to keep secondary alignments") and applies no
// mate-consistency filter, because a read whose mate is unmapped can still carry
// a soft-clipped SL. Either filter silently removes part of the signal this
// subworkflow exists to find.
//
// Verified empirically: an identical alignment configuration is not reachable
// with STAR. Across four scoring regimes (--outFilterMismatchNmax 999 / 10 / 10
// with NoverLmax 0.05 / 6 with NoverLmax 0.03) STAR recovered 5-20 SL-bearing
// tails per million reads from a library that contains ~4,000 SL-anchored 5'
// ends, because STAR has no soft-clip penalty parameter and will always absorb
// an 11 nt leader remnant as mismatches instead of clipping it. STAR cannot be
// tuned into this behaviour.
//
// Scoring is supplied by the caller via ext.args so that the scale factor stays
// in run configuration (see workflows/transcriptomics_operons.config):
//     ext.args = '--sp 1,0 --mp 3,1 -k 5 --score-min L,5,-0.4'
//
process HISAT2_ALIGN_SOFTCLIP {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92f546b522b78597047a54609eefcf35a014de082df0e500bc2adf64d2733669/data' :
        'community.wave.seqera.io/library/hisat2_samtools:a0c9b8ccf8116a89' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam")         , emit: bam
    tuple val(meta), path("*.bam.bai")     , emit: bai
    tuple val(meta), path("*.hisat2.log")  , emit: log
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version | sed -n '1s/.*version //p'"), topic: versions, emit: versions_hisat2
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rg     = args.contains("--rg-id") ? "" : "--rg-id ${prefix} --rg SM:${prefix}"
    // meta.single_end ? single reads : interleaved pair of files
    def mate_args = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"

    """
    INDEX=`find -L ./ -name "*.1.ht2*" | sed 's/\\.1.ht2.*\$//'`
    if [ -z "\$INDEX" ]; then
        echo "ERROR: no hisat2 index (*.1.ht2*) found in the task working directory" >&2
        exit 1
    fi

    # HISAT2 creates named FIFOs in /tmp for its threaded gzip input, naming them
    # from its own PID. Under `singularity --pid` every container sees the same
    # low PIDs while sharing the host /tmp, so two concurrent tasks on a node can
    # collide on /tmp/<pid>.inpipe1, and hisat2 dies with
    #     (ERR): mkfifo(/tmp/18.inpipe1) failed.
    # A per-task TMPDIR removes the shared state (samtools sort honours TMPDIR
    # too); --containall in the run config additionally guarantees a private /tmp
    # for anything that ignores TMPDIR.
    export TMPDIR="\${PWD}/tmp"
    mkdir -p "\$TMPDIR"

    hisat2 \\
        -x \$INDEX \\
        ${mate_args} \\
        --threads ${task.cpus} \\
        --temp-directory "\$TMPDIR" \
        --summary-file ${prefix}.hisat2.log \\
        ${rg} \\
        ${args} \\
        2>> ${prefix}.hisat2.log \\
    | samtools sort -m 1G -@ ${task.cpus} -O bam -o ${prefix}.bam -

    samtools index -@ ${task.cpus} ${prefix}.bam
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    echo "" > ${prefix}.hisat2.log
    """
}
