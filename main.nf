nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// genome_annotation.nf — Master Orchestrator for End-to-End Genome Annotation
//
// Complete End-to-End Pipeline:
// 1. REPEAT_ANNOTATION:
//    - Lewis Stevens Centroid Consensus Satellite Profiling (TRF + VSEARCH + SATELLITEPROFILER)
//    - De novo Transposable Element annotation (Earl Grey)
//    - Emits soft-masked assembly for downstream stages
// 2. NCRNA_ANNOTATION:
//    - Infernal cmscan against Rfam covariance models
//    - Generation of structural RNA exclusion mask
// 3. RNASEQ_ALIGN (Realignment on curated assembly):
//    - HISAT2 de novo index building
//    - Realignment of all staged RNA-seq libraries (CRAM/FASTQ) with soft-clip scoring
//    - Unified coordinate-sorted spliced BAM creation
// 4. SPLICED_LEADER_DISCOVERY:
//    - Empirical trans-splicing leader discovery & TSS hint export
// 5. GENE_PREDICTION:
//    - BRAKER3 (realigned RNA-seq BAM + OrthoDB/reference proteins)
//    - Helixer (Deep learning ab initio)
//    - TSEBRA consensus reconciliation & filtering
//    - Protein set BUSCO quality control
// 6. SPLICED_LEADER_OPERONS:
//    - Demarcates polycistronic operons directly against the resulting gene models
// ---------------------------------------------------------------------------

include { GUNZIP                    } from './modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX            } from './modules/nf-core/samtools/faidx/main'
include { REPEAT_ANNOTATION         } from './subworkflows/sanger-tol/repeat_annotation/main'
include { NCRNA_ANNOTATION          } from './subworkflows/sanger-tol/ncrna_annotation/main'
include { RNASEQ_ALIGN              } from './subworkflows/sanger-tol/rnaseq_align/main'
include { SPLICED_LEADER_DISCOVERY  } from './subworkflows/sanger-tol/spliced_leader_discovery/main'
include { GENE_PREDICTION           } from './subworkflows/sanger-tol/gene_prediction/main'
include { FUNCTIONAL_ANNOTATION     } from './subworkflows/sanger-tol/functional_annotation/main'
include { SLOPPR_OPERONS            } from './modules/sanger-tol/sloppr/operons/main'
include { TESCORER                  } from './modules/sanger-tol/tescorer/main'
include { SUBREAD_FEATURECOUNTS     } from './modules/nf-core/subread/featurecounts/main'
include { SAMTOOLS_VIEW as BAM_TO_CRAM } from './modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX as INDEX_CRAM } from './modules/nf-core/samtools/index/main'

// Top-level Parameters
params.assembly            = params.assembly            ?: null
params.fai                 = params.fai                 ?: null
params.outdir              = params.outdir              ?: "${launchDir}/results_annotation"
params.interproscan_dir    = params.interproscan_dir    ?: null
params.skip_functional_annotation = params.skip_functional_annotation ?: false
params.skip_te_scoring     = params.skip_te_scoring     ?: false
params.orthogroups_tsv     = params.orthogroups_tsv     ?: null
params.satellites_bed      = params.satellites_bed      ?: null

// Repeat parameters
params.window_size         = params.window_size         ?: 500000
params.min_period          = params.min_period          ?: 50
params.min_copies          = params.min_copies          ?: 100
params.min_identity        = params.min_identity        ?: 85
params.clade               = params.clade               ?: 'rhabditida'
params.skip_repeats        = params.skip_repeats        ?: false
params.skip_te             = params.skip_te             ?: false
params.earlgrey_dir        = params.earlgrey_dir        ?: null
params.earlgrey_seed_dir   = params.earlgrey_seed_dir   ?: null

// ncRNA parameters
params.rfam_cm             = params.rfam_cm             ?: null
params.ncrna_gff           = params.ncrna_gff           ?: null
params.skip_ncrna          = params.skip_ncrna          ?: false

// RNA-seq Realignment & Samplesheet parameters
params.rnaseq_samplesheet  = params.rnaseq_samplesheet  ?: null
params.bam                 = params.bam                 ?: null
params.bams                = params.bams                ?: null

// Spliced Leader parameters
params.skip_sl_discovery   = params.skip_sl_discovery   ?: false
params.min_support         = params.min_support         ?: 5
params.max_dist            = params.max_dist            ?: 1000

// Gene Prediction & Pre-computed Ingestion parameters
params.proteins            = params.proteins            ?: null
params.helixer_lineage     = params.helixer_lineage     ?: 'invertebrate'
params.busco_lineage       = params.busco_lineage       ?: 'chromadorea_odb12'
params.busco_lineages_path = params.busco_lineages_path ?: null
params.skip_gene_prediction = params.skip_gene_prediction ?: false

// Pre-computed bypass parameters
params.braker_gtf          = params.braker_gtf          ?: null
params.braker_aa           = params.braker_aa           ?: null
params.braker_hints        = params.braker_hints        ?: null
params.helixer_gff3        = params.helixer_gff3        ?: null
params.helixer_gtf         = params.helixer_gtf         ?: null
params.interproscan_tsv    = params.interproscan_tsv    ?: null
params.interproscan_gff3   = params.interproscan_gff3   ?: null

workflow {
    def commit_suffix = workflow.commitId ? workflow.commitId.take(9) : 'dev'

    if (!params.assembly) {
        error "Mandatory parameter --assembly <path/to/assembly.fa> not specified."
    }

    ch_assembly_file = file(params.assembly, checkIfExists: true)
    ch_assembly_meta = [ id: ch_assembly_file.simpleName ]

    // Transparent on-the-fly decompression for gzipped assemblies (.fa.gz, .fasta.gz)
    if (params.assembly.toString().endsWith('.gz')) {
        GUNZIP(
            channel.of([ ch_assembly_meta, ch_assembly_file ])
        )
        ch_assembly = GUNZIP.out.gunzip
    } else {
        ch_assembly = channel.of([ ch_assembly_meta, ch_assembly_file ])
    }

    // FASTA Index
    if (params.fai && !params.assembly.toString().endsWith('.gz')) {
        ch_fai = channel.of(file(params.fai, checkIfExists: true))
    } else {
        SAMTOOLS_FAIDX(
            ch_assembly.map { meta, fasta -> [ meta, fasta, [] ] },
            false
        )
        ch_fai = SAMTOOLS_FAIDX.out.fai.map { _meta, fai -> fai }
    }

    // -----------------------------------------------------------------------
    // Stage 1: Repeats & Satellites (Lewis Stevens Centroids + Earl Grey)
    // -----------------------------------------------------------------------
    ch_masked_assembly = ch_assembly
    ch_sat_masked_assembly = ch_assembly

    if (!params.skip_repeats) {
        def effective_eg_dir = params.earlgrey_seed_dir ?: params.earlgrey_dir
        REPEAT_ANNOTATION(
            ch_assembly,
            ch_fai,
            params.window_size,
            params.min_period,
            params.min_copies,
            params.min_identity,
            params.clade,
            params.skip_te,
            effective_eg_dir
        )
        ch_masked_assembly     = REPEAT_ANNOTATION.out.masked_fasta
        ch_sat_masked_assembly = REPEAT_ANNOTATION.out.sat_masked_fasta
    }

    // -----------------------------------------------------------------------
    // Stage 2: ncRNA Structural Annotation (Infernal cmscan)
    // -----------------------------------------------------------------------
    ch_ncrna_gff = channel.empty()

    if (params.ncrna_gff) {
        ch_ncrna_gff = channel.value(file(params.ncrna_gff, checkIfExists: true))
    } else if (!params.skip_ncrna && params.rfam_cm) {
        ch_cmfile = file(params.rfam_cm, checkIfExists: true)
        NCRNA_ANNOTATION(
            ch_assembly,
            ch_cmfile
        )
        ch_ncrna_gff = NCRNA_ANNOTATION.out.gff3.map { _meta, gff -> gff }
    }

    // -----------------------------------------------------------------------
    // Stage 3: RNA-seq Realignment (HISAT2 de novo index + Soft-clip alignment)
    // -----------------------------------------------------------------------
    ch_primary_bam = channel.empty()
    ch_individual_bams = channel.empty()

    if (params.rnaseq_samplesheet) {
        ch_samplesheet = channel.fromPath(params.rnaseq_samplesheet, checkIfExists: true)
        RNASEQ_ALIGN(
            ch_assembly,
            ch_fai,
            ch_samplesheet
        )
        ch_primary_bam = RNASEQ_ALIGN.out.merged_bam.map { _meta, bam -> bam }
        ch_individual_bams = RNASEQ_ALIGN.out.individual_bams

        // Convert merged BAM to reference-aligned CRAM and index
        ch_bam_to_cram = RNASEQ_ALIGN.out.merged_bam.map { meta, bam -> [ meta, bam, [] ] }
        ch_fasta_fai   = ch_assembly.join(ch_fai.map { fai -> [ [ id: ch_assembly_meta.id ], fai ] })
        BAM_TO_CRAM(
            ch_bam_to_cram,
            ch_fasta_fai,
            [ [], [] ],
            [ [], [] ],
            []
        )
        INDEX_CRAM(
            BAM_TO_CRAM.out.cram
        )
    } else if (params.bam) {
        ch_primary_bam = channel.value(file(params.bam, checkIfExists: true))
    }

    // -----------------------------------------------------------------------
    // Stage 4: Spliced Leader Discovery & TSS Hints
    // -----------------------------------------------------------------------
    ch_sl_bams = ch_individual_bams

    if (!params.skip_sl_discovery) {
        if (!params.rnaseq_samplesheet && (params.bams || params.bam)) {
            def raw_bams_arg = params.bams ?: params.bam
            def bam_file_list = raw_bams_arg.toString().contains(',')
                ? raw_bams_arg.tokenize(',').collect { fstr -> file(fstr.trim(), checkIfExists: true) }
                : file(raw_bams_arg, checkIfExists: true)

            ch_sl_bams = (bam_file_list instanceof List ? channel.fromList(bam_file_list) : channel.fromPath(raw_bams_arg, checkIfExists: true))
                .map { bfile ->
                    def meta = [ id: bfile.baseName.replaceAll('\\.slidr$', '').replaceAll('\\.sorted$', '') ]
                    def bai = file("${bfile}.bai")
                    if (!bai.exists()) {
                        bai = file("${bfile.toString().replaceAll('\\.bam$', '')}.bai")
                    }
                    [ meta, bfile, bai ]
                }
        }

        SPLICED_LEADER_DISCOVERY(
            ch_sl_bams,
            ch_assembly,
            ch_ncrna_gff,
            params.min_support,
            10, 'AG', params.max_dist, 'GT',
            '.{10,60}(AATTTTTGGA|AATTTTGGA|AATTTGGA|AT{4,6}G)',
            120, 3, 10, 3, 1, 500
        )
    }

    // -----------------------------------------------------------------------
    // Stage 5: Dual Gene Prediction (BRAKER3 + Helixer + TSEBRA consensus)
    // -----------------------------------------------------------------------
    ch_consensus_gtf = channel.empty()

    if (!params.skip_gene_prediction) {
        ch_proteins       = params.proteins ? file(params.proteins, checkIfExists: true) : []
        ch_busco_lineages = params.busco_lineages_path ? file(params.busco_lineages_path, checkIfExists: true) : []

        file_braker_gtf   = params.braker_gtf   ? file(params.braker_gtf, checkIfExists: true)   : null
        file_braker_aa    = params.braker_aa    ? file(params.braker_aa, checkIfExists: true)    : null
        file_braker_hints = params.braker_hints ? file(params.braker_hints, checkIfExists: true) : null
        file_helixer_gff3 = params.helixer_gff3 ? file(params.helixer_gff3, checkIfExists: true) : null
        file_helixer_gtf  = params.helixer_gtf  ? file(params.helixer_gtf, checkIfExists: true)  : null

        GENE_PREDICTION(
            ch_sat_masked_assembly,
            ch_primary_bam,
            ch_proteins,
            ch_busco_lineages,
            params.helixer_lineage,
            params.busco_lineage,
            file_braker_gtf,
            file_braker_aa,
            file_braker_hints,
            file_helixer_gff3,
            file_helixer_gtf
        )
        ch_consensus_gtf = GENE_PREDICTION.out.consensus_gtf

        // -------------------------------------------------------------------
        // Stage 6: Functional Annotation (InterProScan)
        // -------------------------------------------------------------------
        if (!params.skip_functional_annotation) {
            if (params.interproscan_tsv) {
                // Pre-computed InterProScan results provided: publish/expose directly
                file_ipr_tsv  = file(params.interproscan_tsv, checkIfExists: true)
                file_ipr_gff3 = params.interproscan_gff3 ? file(params.interproscan_gff3, checkIfExists: true) : null
                ch_assembly.map { meta, _f ->
                    def outdir_ipr = file("${params.outdir}/functional_annotation/interproscan")
                    outdir_ipr.mkdirs()
                    // Create convenient symlinks/copies in standard output location
                    def target_tsv = file("${outdir_ipr}/${meta.id}.iprscan.tsv")
                    if (!target_tsv.exists()) {
                        file_ipr_tsv.copyTo(target_tsv)
                    }
                    if (file_ipr_gff3) {
                        def target_gff3 = file("${outdir_ipr}/${meta.id}.iprscan.gff3")
                        if (!target_gff3.exists()) {
                            file_ipr_gff3.copyTo(target_gff3)
                        }
                    }
                }
            } else if (params.interproscan_dir) {
                ch_ipr_dir = file(params.interproscan_dir, checkIfExists: true)
                ch_proteins_to_annotate = GENE_PREDICTION.out.representative_aa

                FUNCTIONAL_ANNOTATION(
                    ch_proteins_to_annotate,
                    ch_ipr_dir
                )
            }
        }
    }

    // -----------------------------------------------------------------------
    // Stage 7: Operon Demarcation (SLOPPR)
    // -----------------------------------------------------------------------
    ch_operons_gff3 = channel.empty()
    if (!params.skip_sl_discovery && !params.skip_gene_prediction) {
        ch_representative_gtf = GENE_PREDICTION.out.representative_gtf.map { _meta, gtf -> gtf }.first()

        SLOPPR_OPERONS(
            SPLICED_LEADER_DISCOVERY.out.sl_bed,
            ch_representative_gtf,
            params.max_dist
        )
        ch_operons_gff3 = SLOPPR_OPERONS.out.operons_gff
    }

    // -----------------------------------------------------------------------
    // Stage 8: Transposable Element (TE) Multi-Evidence Scoring
    // -----------------------------------------------------------------------
    if (!params.skip_gene_prediction && !params.skip_te_scoring) {
        ch_te_gff3 = (!params.skip_repeats) ? REPEAT_ANNOTATION.out.te_gff3 : channel.empty()

        ch_sat_bed = params.satellites_bed ?
            channel.of([ [ id: ch_assembly_meta.id ], file(params.satellites_bed, checkIfExists: true) ]) :
            ((!params.skip_repeats) ? REPEAT_ANNOTATION.out.satellite_bed : channel.empty())

        ch_orthogroups = params.orthogroups_tsv ?
            channel.of([ [ id: ch_assembly_meta.id ], file(params.orthogroups_tsv, checkIfExists: true) ]) :
            channel.empty()

        ch_ipr_tsv = params.interproscan_tsv ?
            channel.of([ [ id: ch_assembly_meta.id ], file(params.interproscan_tsv, checkIfExists: true) ]) :
            (params.interproscan_dir ? FUNCTIONAL_ANNOTATION.out.tsv : channel.empty())

        // Combine all evidence tracks for the representative models
        // Remap all supplementary tracks to match the primary assembly meta ID
        ch_scorer_input = GENE_PREDICTION.out.representative_gtf
            .join(GENE_PREDICTION.out.representative_aa)
            .join(ch_ipr_tsv.map { _m, f -> [ [ id: ch_assembly_meta.id ], f ] }, remainder: true)
            .join(ch_te_gff3.map { _m, f -> [ [ id: ch_assembly_meta.id ], f ] }, remainder: true)
            .join(ch_sat_bed.map { _m, f -> [ [ id: ch_assembly_meta.id ], f ] }, remainder: true)
            .join(ch_operons_gff3.map { _m, f -> [ [ id: ch_assembly_meta.id ], f ] }, remainder: true)
            .join(ch_orthogroups.map { _m, f -> [ [ id: ch_assembly_meta.id ], f ] }, remainder: true)
            .filter { item ->
                // Guard: Only process the record whose meta corresponds to the assembly and has valid gtf/aa
                item[0]?.id == ch_assembly_meta.id && item[1] != null && item[2] != null
            }
            .map { item ->
                def meta = item[0]
                def gtf  = item[1]
                def aa   = item[2]
                def ipr  = item.size() > 3 ? item[3] : []
                def eg   = item.size() > 4 ? item[4] : []
                def sat  = item.size() > 5 ? item[5] : []
                def op   = item.size() > 6 ? item[6] : []
                def og   = item.size() > 7 ? item[7] : []
                [ meta, gtf, aa, ipr ?: [], eg ?: [], sat ?: [], op ?: [], og ?: [] ]
            }

        TESCORER(ch_scorer_input)
    }

    // -----------------------------------------------------------------------
    // Stage 9: RNA-seq Gene Abundance Estimation (featureCounts)
    // -----------------------------------------------------------------------
    if (params.rnaseq_samplesheet && !params.skip_gene_prediction) {
        ch_counts_input = RNASEQ_ALIGN.out.merged_bam
            .map { meta, bam -> [ [ id: ch_assembly_meta.id, single_end: false ], bam ] }
            .join(GENE_PREDICTION.out.representative_gtf)
            .map { meta, bam, gtf -> [ meta, [ bam ], gtf ] }

        SUBREAD_FEATURECOUNTS(
            ch_counts_input
        )
    }
}
