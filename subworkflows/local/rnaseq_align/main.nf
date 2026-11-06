//
// Subworkflow: RNASEQ_ALIGN
// Build HISAT2 index, convert unaligned CRAMs to FASTQs, align with soft-clip permissive scoring,
// convert/merge pre-aligned CRAMs and BAMs, and merge into a unified coordinate-sorted BAM.
//

include { HISAT2_BUILD                                                    } from '../../../modules/nf-core/hisat2/build/main'
include { HISAT2_ALIGN_SOFTCLIP as SLIDR_ALIGNREADS                       } from '../../../modules/local/slidr/alignreads/main'
include { SAMTOOLS_INDEX                                                  } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                                                  } from '../../../modules/local/samtools/merge/main'
include { SAMTOOLS_COLLATEFASTQ                                           } from '../../../modules/nf-core/samtools/collatefastq/main'

def resolve_path(p) {
    if (!p) return null
    def f = file(p)
    if (f.exists()) return f
    if (params.rnaseq_samplesheet) {
        def rel = file("${file(params.rnaseq_samplesheet).parent}/${p}")
        if (rel.exists()) return rel
    }
    return file(p, checkIfExists: true)
}

def is_prealigned_cram(row) {
    return (row.aligned_cram != null && row.aligned_cram != '') ||
           (row.cram && (row.aligned?.toString()?.toLowerCase() == 'true' || row.is_aligned?.toString()?.toLowerCase() == 'true'))
}

process CRAM_TO_BAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/../../../modules/nf-core/samtools/view/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c2f02ceb8cbe5aa37930dcbc0de6f1d07dfb292178c7e780c44ab2f556bf829f/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--8a3ddf262165ba23'}"

    input:
    tuple val(meta), path(cram), path(crai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam_bai
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view \\
        -@ ${task.cpus} \\
        -b \\
        -T "${fasta}" \\
        -o "${prefix}.bam##idx##${prefix}.bam.bai" \\
        --write-index \\
        "${cram}"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    """
}

workflow RNASEQ_ALIGN {
    take:
    ch_assembly       // channel: [ val(meta), path(fasta) ]
    ch_fai            // channel: path(fai)
    ch_samplesheet    // channel: path(samplesheet_csv)

    main:

    // 1. Build HISAT2 index de novo
    HISAT2_BUILD(
        ch_assembly.map { meta, fasta -> [ meta, fasta, [], [] ] },
        false
    )
    ch_hisat2_index = HISAT2_BUILD.out.index.first()

    // 2. Parse Samplesheet
    ch_rows = ch_samplesheet
        .splitCsv(header: true)

    // Pre-aligned CRAM rows: converted directly without realignment
    ch_cram_prealigned = ch_rows
        .filter { row -> is_prealigned_cram(row) }
        .map { row ->
            def single_end = (row.single_end != null && row.single_end != '') ? row.single_end.toBoolean() : false
            def strandedness = row.strandedness ?: 'unstranded'
            def meta = [ id: row.sample, single_end: single_end, strandedness: strandedness ]
            def cram = row.aligned_cram ? resolve_path(row.aligned_cram) : resolve_path(row.cram)
            def crai = row.aligned_crai ? resolve_path(row.aligned_crai) :
                       (row.crai ? resolve_path(row.crai) :
                       (file("${cram}.crai").exists() ? file("${cram}.crai") : []))
            [ meta, cram, crai ]
        }

    // Raw unaligned CRAMs: collate to FASTQ then align with HISAT2
    ch_cram_input = ch_rows
        .filter { row -> row.cram && !is_prealigned_cram(row) && !row.fastq_1 && !row.bam }
        .map { row ->
            def single_end = (row.single_end != null && row.single_end != '') ? row.single_end.toBoolean() : false
            def strandedness = row.strandedness ?: 'unstranded'
            def meta = [ id: row.sample, single_end: single_end, strandedness: strandedness ]
            [ meta, resolve_path(row.cram) ]
        }

    // FASTQ direct rows
    ch_fastq_direct = ch_rows
        .filter { row -> row.fastq_1 && !row.bam && !row.aligned_cram }
        .map { row ->
            def single_end = (row.single_end != null && row.single_end != '') ? row.single_end.toBoolean() : (row.fastq_2 ? false : true)
            def strandedness = row.strandedness ?: 'unstranded'
            def meta = [ id: row.sample, single_end: single_end, strandedness: strandedness ]
            def reads = single_end ? [ resolve_path(row.fastq_1) ] :
                                     [ resolve_path(row.fastq_1), resolve_path(row.fastq_2) ]
            [ meta, reads ]
        }

    // BAM direct rows
    ch_bam_direct = ch_rows
        .filter { row -> row.bam && !row.aligned_cram }
        .map { row ->
            def single_end = (row.single_end != null && row.single_end != '') ? row.single_end.toBoolean() : false
            def strandedness = row.strandedness ?: 'unstranded'
            def meta = [ id: row.sample, single_end: single_end, strandedness: strandedness ]
            def bam = resolve_path(row.bam)
            def bai = row.bai ? resolve_path(row.bai) : (file("${bam}.bai").exists() ? file("${bam}.bai") : resolve_path("${row.bam}.bai"))
            [ meta, bam, bai ]
        }

    // Convert raw unaligned CRAMs to FASTQ
    ch_empty_ref = channel.of([ [ id: 'empty' ], [], [] ])
    SAMTOOLS_COLLATEFASTQ(
        ch_cram_input,
        ch_empty_ref.first(),
        false
    )

    ch_fastq_from_cram = SAMTOOLS_COLLATEFASTQ.out.fastq.map { meta, fastqs -> [ meta, fastqs ] }
    ch_all_fastqs = ch_fastq_direct.mix(ch_fastq_from_cram)

    // 3. Align raw reads using SLIDR soft-clip parameters
    SLIDR_ALIGNREADS(
        ch_all_fastqs,
        ch_hisat2_index
    )

    // Convert pre-aligned CRAMs directly to BAM + BAI
    ch_assembly_fasta = ch_assembly.map { _meta, fasta -> fasta }.first()
    ch_assembly_fai   = ch_fai.first()

    CRAM_TO_BAM(
        ch_cram_prealigned,
        ch_assembly_fasta,
        ch_assembly_fai
    )

    ch_aligned_bams = SLIDR_ALIGNREADS.out.bam
        .join(SLIDR_ALIGNREADS.out.bai)
        .mix(ch_bam_direct)
        .mix(CRAM_TO_BAM.out.bam_bai)

    // 4. Merge all BAMs into unified BAM for BRAKER3
    ch_all_bam_files = ch_aligned_bams.map { _meta, bam, _bai -> bam }.collect()
    ch_all_bai_files = ch_aligned_bams.map { _meta, _bam, bai -> bai }.collect()

    ch_merge_meta = [ id: 'rnaseq_merged_spliced' ]
    ch_merge_input = ch_all_bam_files
        .toList()
        .combine(ch_all_bai_files.toList())
        .map { bams, bais -> [ ch_merge_meta, bams, bais ] }

    SAMTOOLS_MERGE(
        ch_merge_input,
        ch_assembly.combine(ch_fai).map { meta, fasta, fai -> [ meta, fasta, fai, [] ] }.first(),
        'bai'
    )

    emit:
    merged_bam      = SAMTOOLS_MERGE.out.bam          // [ meta, rnaseq_merged_spliced.bam ]
    merged_bai      = SAMTOOLS_MERGE.out.index        // [ meta, rnaseq_merged_spliced.bam.bai ]
    individual_bams = ch_aligned_bams                 // [ meta, bam, bai ]
}
