/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW: NEMATODEANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Dual-evidence protein + RNA tournament annotation pipeline tailored for nematodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Subworkflows
include { REPEAT_ANNOTATION         } from '../subworkflows/local/repeat_annotation/main'
include { RNASEQ_ALIGN              } from '../subworkflows/local/rnaseq_align/main'
include { BAM_TO_CRAM               } from '../subworkflows/local/bam_to_cram/main'
include { SPLICED_LEADER_DISCOVERY  } from '../subworkflows/local/spliced_leader_discovery/main'
include { GENE_PREDICTION           } from '../subworkflows/local/gene_prediction/main'
include { FUNCTIONAL_ANNOTATION     } from '../subworkflows/local/functional_annotation/main'

// Modules
include { SUBREAD_FEATURECOUNTS     } from '../modules/nf-core/subread/featurecounts/main'

workflow NEMATODEANNOTATION {
    take:
    ch_assembly            // channel: [ val(meta), path(fasta) ]
    ch_fai                 // channel: [ val(meta), path(fai) ] or path(fai)
    ch_rnaseq_samplesheet  // path(samplesheet) or null
    ch_proteins            // path(proteins_fasta) or []
    ch_busco_lineages      // path(busco_lineages_dir) or []

    main:
    ch_versions = channel.empty()

    // 1. Repeat Annotation (Centroid Consensus Satellites + Earl Grey TEs)
    ch_fai_plain = ch_fai.map { entry -> entry instanceof List ? (entry.size() > 1 ? entry[1] : entry[0]) : entry }
    def skip_te_effective = params.skip_te || params.skip_earlgrey

    REPEAT_ANNOTATION(
        ch_assembly,
        ch_fai_plain,
        params.window_size ?: 500000,
        params.min_period ?: 50,
        params.min_copies ?: 100,
        params.min_identity ?: 85,
        params.clade ?: 'rhabditida',
        skip_te_effective,
        params.earlgrey_dir ?: null,
        params.earlgrey_seed_dir ?: null,
        channel.empty()
    )
    ch_versions = ch_versions.mix(REPEAT_ANNOTATION.out.versions)

    def ch_assembly_for_genes = (params.use_sat_only_mask_for_genes || skip_te_effective) ?
        REPEAT_ANNOTATION.out.sat_masked_fasta : REPEAT_ANNOTATION.out.masked_fasta

    // 2. RNA-seq Alignment via HISAT2 & Reference CRAM Conversion
    ch_cram = channel.empty()
    ch_crai = channel.empty()
    ch_bam_for_gene_pred = channel.empty()

    if (ch_rnaseq_samplesheet) {
        ch_samplesheet = channel.fromPath(ch_rnaseq_samplesheet, checkIfExists: true)

        RNASEQ_ALIGN(
            ch_assembly,
            ch_fai_plain,
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(RNASEQ_ALIGN.out.versions)

        // 3. Trans-spliced Leader (SL) Tail Extraction & Profiling
        SPLICED_LEADER_DISCOVERY(
            RNASEQ_ALIGN.out.individual_bams,
            ch_assembly,
            channel.empty(),
            params.min_support ?: 5,
            10, 'AG', params.max_dist ?: 1000, 'GT',
            '.{10,60}(AATTTTTGGA|AATTTTGGA|AATTTGGA|AT{4,6}G)',
            120, 3, 10, 3, 1, 500
        )
        ch_versions = ch_versions.mix(SPLICED_LEADER_DISCOVERY.out.versions)

        // 4. Convert Merged RNA-seq BAM to Assembly-Referenced CRAM
        BAM_TO_CRAM(
            RNASEQ_ALIGN.out.merged_bam,
            ch_assembly,
            ch_fai_plain
        )
        ch_cram = BAM_TO_CRAM.out.cram
        ch_crai = BAM_TO_CRAM.out.crai
        ch_versions = ch_versions.mix(BAM_TO_CRAM.out.versions)

        // Pass [bam, bai] tuple to GENE_PREDICTION
        ch_bam_for_gene_pred = RNASEQ_ALIGN.out.merged_bam
            .join(RNASEQ_ALIGN.out.merged_bai)
            .map { _meta, bam, bai -> [ bam, bai ] }
    }

    // 5. Dual-Evidence Gene Prediction Multi-Candidate Tournament
    def skip_helixer_val = params.skip_helixer || !params.run_helixer

    GENE_PREDICTION(
        ch_assembly_for_genes,
        ch_bam_for_gene_pred,
        ch_proteins,
        ch_busco_lineages,
        params.helixer_lineage ?: 'invertebrate',
        params.busco_lineage   ?: 'chromadorea_odb12',
        params.braker_gtf      ? file(params.braker_gtf, checkIfExists: true) : null,
        params.braker_aa       ? file(params.braker_aa, checkIfExists: true) : null,
        params.braker_hints    ? file(params.braker_hints, checkIfExists: true) : null,
        params.helixer_gff3    ? file(params.helixer_gff3, checkIfExists: true) : null,
        params.helixer_gtf     ? file(params.helixer_gtf, checkIfExists: true) : null,
        skip_helixer_val
    )
    ch_versions = ch_versions.mix(GENE_PREDICTION.out.versions)

    // 6. Functional Annotation of Winning Model via InterProScan
    if (!params.skip_functional_annotation && params.interproscan_dir) {
        FUNCTIONAL_ANNOTATION(
            GENE_PREDICTION.out.representative_faa,
            file(params.interproscan_dir, checkIfExists: true)
        )
        ch_versions = ch_versions.mix(FUNCTIONAL_ANNOTATION.out.versions)
    }

    // 7. RNA-seq Gene Abundance Estimation via featureCounts
    if (ch_rnaseq_samplesheet) {
        ch_counts_input = RNASEQ_ALIGN.out.merged_bam
            .combine(GENE_PREDICTION.out.representative_gtf)
            .map { aln_meta, bam, _gtf_meta, gtf ->
                def counts_meta = [
                    id: aln_meta.id,
                    single_end: aln_meta.single_end ?: false,
                    strandedness: aln_meta.strandedness ?: 'unstranded'
                ]
                [ counts_meta, [ bam ], gtf ]
            }

        SUBREAD_FEATURECOUNTS(
            ch_counts_input
        )
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions_subread)
    }

    emit:
    representative_gtf              = GENE_PREDICTION.out.representative_gtf
    representative_faa              = GENE_PREDICTION.out.representative_faa
    representative_fna              = GENE_PREDICTION.out.representative_fna
    decision_tsv                    = GENE_PREDICTION.out.representative_tsv
    representative_busco_full_table = GENE_PREDICTION.out.representative_busco_full_table
    suboptimal_models               = GENE_PREDICTION.out.suboptimal_models_archive
    cram                            = ch_cram
    crai                            = ch_crai
    versions                        = ch_versions
}
