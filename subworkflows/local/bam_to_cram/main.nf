//
// Subworkflow: BAM_TO_CRAM
// Convert coordinate-sorted BAM to reference-aligned CRAM and index
//

include { SAMTOOLS_INDEX as INDEX_CRAM } from '../../../modules/nf-core/samtools/index/main'

process BAM_TO_CRAM_PROCESS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/../../../modules/nf-core/samtools/view/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c2f02ceb8cbe5aa37930dcbc0de6f1d07dfb292178c7e780c44ab2f556bf829f/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--8a3ddf262165ba23'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${prefix}.cram"), emit: cram
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--output-fmt cram'
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view \\
        -@ ${task.cpus} \\
        ${args} \\
        --reference "${fasta}" \\
        -o "${prefix}.cram" \\
        "${bam}"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cram
    """
}

workflow BAM_TO_CRAM {
    take:
    ch_bam        // channel: [ val(meta), path(bam), path(bai) ] or [ val(meta), path(bam) ]
    ch_fasta      // channel: [ val(meta_fasta), path(fasta) ] or path(fasta)
    ch_fai        // channel: [ val(meta_fai), path(fai) ] or path(fai)

    main:
    ch_versions = channel.empty()

    // Normalize BAM input channel to [ meta, bam, bai ]
    ch_bam_normalized = ch_bam.map { entry ->
        if (entry instanceof List) {
            if (entry.size() == 2) {
                def (meta, bam) = entry
                return [ meta, bam, [] ]
            } else if (entry.size() >= 3) {
                def (meta, bam, bai) = entry
                return [ meta, bam, bai ?: [] ]
            }
        }
        return [ [ id: file(entry).baseName ], entry, [] ]
    }

    // Normalize reference FASTA and FAI to file paths
    ch_fasta_file = ch_fasta.map { entry ->
        entry instanceof List ? (entry.size() > 1 ? entry[1] : entry[0]) : entry
    }.first()

    ch_fai_file = ch_fai.map { entry ->
        entry instanceof List ? (entry.size() > 1 ? entry[1] : entry[0]) : entry
    }.first()

    // 1. Convert BAM to CRAM with reference FASTA
    BAM_TO_CRAM_PROCESS(
        ch_bam_normalized,
        ch_fasta_file,
        ch_fai_file
    )
    ch_versions = ch_versions.mix(BAM_TO_CRAM_PROCESS.out.versions_samtools)

    // 2. Index CRAM to generate .crai
    INDEX_CRAM(
        BAM_TO_CRAM_PROCESS.out.cram
    )
    ch_versions = ch_versions.mix(INDEX_CRAM.out.versions_samtools)

    ch_cram_crai = BAM_TO_CRAM_PROCESS.out.cram
        .join(INDEX_CRAM.out.index)
        .map { meta, cram, crai -> [ meta, cram, crai ] }

    emit:
    cram      = BAM_TO_CRAM_PROCESS.out.cram    // [ val(meta), path(cram) ]
    crai      = INDEX_CRAM.out.index            // [ val(meta), path(crai) ]
    cram_crai = ch_cram_crai                    // [ val(meta), path(cram), path(crai) ]
    versions  = ch_versions                     // topic: versions
}
