//
// Subworkflow: RNASEQ_ALIGN
// Build HISAT2 index, convert unaligned CRAMs to FASTQs, align with soft-clip permissive scoring,
// and merge into a unified coordinate-sorted BAM.
//

include { HISAT2_BUILD                                                    } from '../../../modules/nf-core/hisat2/build/main'
include { HISAT2_ALIGN_SOFTCLIP as SLIDR_ALIGNREADS                       } from '../../../modules/local/slidr/alignreads/main'
include { SAMTOOLS_INDEX                                                  } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                                                  } from '../../../modules/local/samtools/merge/main'
include { SAMTOOLS_COLLATEFASTQ                                           } from '../../../modules/nf-core/samtools/collatefastq/main'

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

    ch_cram_input = ch_rows
        .filter { row -> row.cram && !row.fastq_1 && !row.bam }
        .map { row ->
            def meta = [ id: row.sample, single_end: false ]
            [ meta, file(row.cram, checkIfExists: true) ]
        }

    ch_fastq_direct = ch_rows
        .filter { row -> row.fastq_1 && !row.bam }
        .map { row ->
            def single_end = row.fastq_2 ? false : true
            def meta = [ id: row.sample, single_end: single_end ]
            def reads = single_end ? [ file(row.fastq_1, checkIfExists: true) ] :
                                     [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]
            [ meta, reads ]
        }

    ch_bam_direct = ch_rows
        .filter { row -> row.bam }
        .map { row ->
            def meta = [ id: row.sample ]
            def bam = file(row.bam, checkIfExists: true)
            def bai = row.bai ? file(row.bai, checkIfExists: true) : file("${row.bam}.bai", checkIfExists: true)
            [ meta, bam, bai ]
        }

    // Convert CRAMs to FASTQ (pass empty reference tuple so unmapped reads without reference slices are collated directly)
    ch_empty_ref = channel.of([ [ id: 'empty' ], [], [] ])
    SAMTOOLS_COLLATEFASTQ(
        ch_cram_input,
        ch_empty_ref.first(),
        false
    )

    ch_fastq_from_cram = SAMTOOLS_COLLATEFASTQ.out.fastq.map { meta, fastqs -> [ meta, fastqs ] }
    ch_all_fastqs = ch_fastq_direct.mix(ch_fastq_from_cram)

    // 3. Align reads using SLIDR soft-clip parameters
    SLIDR_ALIGNREADS(
        ch_all_fastqs,
        ch_hisat2_index
    )

    ch_aligned_bams = SLIDR_ALIGNREADS.out.bam
        .join(SLIDR_ALIGNREADS.out.bai)
        .mix(ch_bam_direct)

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
