/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW: NEMATODEANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Dual-evidence protein + RNA tournament annotation pipeline tailored for nematodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Subworkflows
include { REPEAT_MASKING_SATELLITES } from '../subworkflows/sanger-tol/repeat_masking_satellites/main'
include { REPEAT_ANNOTATION         } from '../subworkflows/sanger-tol/repeat_annotation/main'
include { RNASEQ_ALIGN              } from '../subworkflows/sanger-tol/rnaseq_align/main'
include { SPLICED_LEADER_DISCOVERY  } from '../subworkflows/sanger-tol/spliced_leader_discovery/main'
include { GENE_PREDICTION           } from '../subworkflows/sanger-tol/gene_prediction/main'
include { FUNCTIONAL_ANNOTATION     } from '../subworkflows/sanger-tol/functional_annotation/main'

// Modules
include { SAMTOOLS_VIEW as BAM_TO_CRAM } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX as INDEX_CRAM } from '../modules/nf-core/samtools/index/main'
include { SUBREAD_FEATURECOUNTS        } from '../modules/nf-core/subread/featurecounts/main'

workflow NEMATODEANNOTATION {
    take:
    ch_assembly            // channel: [ val(meta), path(fasta) ]
    ch_fai                 // channel: [ val(meta), path(fai) ]
    ch_rnaseq_samplesheet  // path(samplesheet) or null
    ch_proteins            // path(proteins_fasta) or []
    ch_busco_lineages      // path(busco_lineages_dir) or []

    main:
    ch_versions = channel.empty()

    // 1. Satellite & Tandem Repeat Hierarchical Masking
    REPEAT_MASKING_SATELLITES(
        ch_assembly
    )
    ch_versions = ch_versions.mix(REPEAT_MASKING_SATELLITES.out.versions)

    // 2. Transposable Element Annotation via Earl Grey
    REPEAT_ANNOTATION(
        REPEAT_MASKING_SATELLITES.out.masked_fasta,
        params.species_name ?: "Nematoda"
    )
    ch_versions = ch_versions.mix(REPEAT_ANNOTATION.out.versions)

    // 3. RNA-seq Alignment via HISAT2
    ch_cram = channel.empty()
    ch_merged_bam = channel.empty()

    if (ch_rnaseq_samplesheet) {
        ch_samplesheet = channel.fromPath(ch_rnaseq_samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def meta = [ id: row.sample, single_end: false, strandedness: row.strandedness ?: 'unstranded' ]
                def cram = file(row.cram, checkIfExists: true)
                def crai = row.crai ? file(row.crai, checkIfExists: true) : file("${row.cram}.crai")
                [ meta, cram, crai ]
            }

        RNASEQ_ALIGN(
            ch_samplesheet,
            ch_assembly,
            ch_fai
        )
        ch_versions = ch_versions.mix(RNASEQ_ALIGN.out.versions)

        // 4. Trans-spliced Leader (SL) Tail Extraction & Profiling
        SPLICED_LEADER_DISCOVERY(
            RNASEQ_ALIGN.out.collate_fastq,
            ch_assembly
        )
        ch_versions = ch_versions.mix(SPLICED_LEADER_DISCOVERY.out.versions)

        ch_merged_bam = RNASEQ_ALIGN.out.merged_bam

        // 5. Convert Merged RNA-seq BAM to Assembly-Referenced CRAM
        ch_bam_for_cram = ch_merged_bam.map { meta, bam -> [ meta, bam, [] ] }
        BAM_TO_CRAM(
            ch_bam_for_cram,
            ch_assembly,
            ch_fai
        )
        ch_versions = ch_versions.mix(BAM_TO_CRAM.out.versions_samtools)

        INDEX_CRAM(
            BAM_TO_CRAM.out.cram
        )
        ch_versions = ch_versions.mix(INDEX_CRAM.out.versions_samtools)
    }

    // 6. Dual-Evidence Gene Prediction Multi-Candidate Tournament
    ch_bam_for_gene_pred = ch_rnaseq_samplesheet ? ch_merged_bam.map { _meta, bam -> bam } : channel.empty()

    GENE_PREDICTION(
        REPEAT_ANNOTATION.out.masked_fasta,
        ch_bam_for_gene_pred,
        ch_proteins,
        ch_busco_lineages,
        params.helixer_lineage ?: 'invertebrate',
        params.busco_lineage   ?: 'chromadorea_odb12',
        params.braker_gtf      ? file(params.braker_gtf, checkIfExists: true) : null,
        params.braker_aa       ? file(params.braker_aa, checkIfExists: true) : null,
        params.braker_hints    ? file(params.braker_hints, checkIfExists: true) : null,
        params.helixer_gff3    ? file(params.helixer_gff3, checkIfExists: true) : null,
        params.helixer_gtf     ? file(params.helixer_gtf, checkIfExists: true) : null
    )
    ch_versions = ch_versions.mix(GENE_PREDICTION.out.versions)

    // 7. Functional Annotation of Winning Model via InterProScan
    if (!params.skip_functional_annotation && params.interproscan_dir) {
        FUNCTIONAL_ANNOTATION(
            GENE_PREDICTION.out.representative_faa,
            file(params.interproscan_dir, checkIfExists: true)
        )
        ch_versions = ch_versions.mix(FUNCTIONAL_ANNOTATION.out.versions)
    }

    // 8. RNA-seq Gene Abundance Estimation via featureCounts
    if (ch_rnaseq_samplesheet) {
        ch_counts_input = ch_merged_bam
            .join(GENE_PREDICTION.out.representative_gtf)
            .map { meta, bam, gtf ->
                [ meta, [ bam ], gtf ]
            }

        SUBREAD_FEATURECOUNTS(
            ch_counts_input
        )
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions_subread)
    }

    emit:
    representative_gtf = GENE_PREDICTION.out.representative_gtf
    representative_faa = GENE_PREDICTION.out.representative_faa
    representative_fna = GENE_PREDICTION.out.representative_fna
    decision_tsv       = GENE_PREDICTION.out.representative_tsv
    suboptimal_models  = GENE_PREDICTION.out.suboptimal_models_archive
    versions           = ch_versions
}
