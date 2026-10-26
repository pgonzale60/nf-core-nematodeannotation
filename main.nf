nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// genome_annotation.nf — Master Orchestrator for End-to-End Genome Annotation
//
// Complete End-to-End Pipeline:
// 1. REPEAT_ANNOTATION:
//    - De novo Centroid Consensus Satellite Profiling (TRF + VSEARCH + SATELLITEPROFILER)
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
include { REPEAT_ANNOTATION         } from './subworkflows/local/repeat_annotation/main'
include { NCRNA_ANNOTATION          } from './subworkflows/local/ncrna_annotation/main'
include { NCRNA_TANDEM_FILTER       } from './modules/local/ncrna_tandem_filter/main'
include { RNASEQ_ALIGN              } from './subworkflows/local/rnaseq_align/main'
include { SPLICED_LEADER_DISCOVERY  } from './subworkflows/local/spliced_leader_discovery/main'
include { GENE_PREDICTION           } from './subworkflows/local/gene_prediction/main'
include { FUNCTIONAL_ANNOTATION     } from './subworkflows/local/functional_annotation/main'
include { SLOPPR_OPERONS            } from './modules/local/sloppr/operons/main'
include { TESCORER                  } from './modules/local/tescorer/main'
include { SUBREAD_FEATURECOUNTS     } from './modules/nf-core/subread/featurecounts/main'
include { SAMTOOLS_VIEW as BAM_TO_CRAM } from './modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX as INDEX_CRAM   } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_PRECOMPUTED_BAM } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_SL_BAM } from './modules/nf-core/samtools/index/main'
include { VALIDATE_RNASEQ_ALIGNMENT      } from './modules/local/validate_rnaseq_alignment/main'
include { GFFREAD_PROTEIN           } from './modules/local/gffread/main'
include { GFFREAD_CDS               } from './modules/local/gffread/main'
include { softwareVersionsToYAML    } from './subworkflows/nf-core/utils_nfcore_pipeline/main'
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

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
params.structural_mask_bed = params.structural_mask_bed ?: null

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
params.skip_helixer        = params.skip_helixer        ?: false
params.run_helixer         = params.run_helixer         ?: false
params.busco_lineage       = params.busco_lineage       ?: 'chromadorea_odb12'
params.busco_lineages_path = params.busco_lineages_path ?: null
params.skip_gene_prediction = params.skip_gene_prediction ?: false

// Pre-computed bypass parameters
params.representative_gtf   = params.representative_gtf   ?: null
params.representative_faa   = params.representative_faa   ?: null
params.braker_gtf          = params.braker_gtf          ?: null
params.braker_aa           = params.braker_aa           ?: null
params.braker_hints        = params.braker_hints        ?: null
params.helixer_gff3        = params.helixer_gff3        ?: null
params.helixer_gtf         = params.helixer_gtf         ?: null
params.interproscan_tsv    = params.interproscan_tsv    ?: null
params.interproscan_gff3   = params.interproscan_gff3   ?: null

workflow {
    if (params.help) {
        log.info paramsHelp()
        exit 0
    }

    // 1. Declarative parameter validation via nf-schema Draft 2020-12
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // 2. Validate RNA-seq samplesheet via nf-schema Draft 2020-12
    if (params.rnaseq_samplesheet) {
        samplesheetToList(params.rnaseq_samplesheet, "${projectDir}/assets/schema_input.json")
    }

    // =========================================================================
    // DOMAIN-SPECIFIC PRE-RUN VERIFICATION GATE
    // (Verifies binary indices & non-empty statuses beyond JSON schema reach)
    // =========================================================================
    def validation_errors = []

    // 1. Assembly non-empty verification
    if (params.assembly) {
        def fa_file = file(params.assembly)
        if (fa_file.isEmpty()) {
            validation_errors.add("Assembly file is empty: ${params.assembly}")
        }
    }

    // 2. ncRNA & Rfam Model Indices Verification
    if (!params.skip_ncrna) {
        if (!params.rfam_cm && !params.ncrna_gff && !params.structural_mask_bed && !params.earlgrey_seed_dir && !params.earlgrey_dir) {
            validation_errors.add("ncRNA annotation is active, but neither --rfam_cm nor --ncrna_gff was provided.")
        } else if (params.rfam_cm) {
            def cm = file(params.rfam_cm)
            if (cm.isEmpty()) {
                validation_errors.add("Rfam covariance model file is empty: ${params.rfam_cm}")
            } else {
                def missing_cmpress = ['.i1m', '.i1i', '.i1f', '.i1p'].findAll { ext -> !file("${params.rfam_cm}${ext}").exists() }
                if (missing_cmpress) {
                    log.warn "[PRE-RUN GATE] Rfam.cm is missing cmpress indices (${missing_cmpress.join(', ')}). Infernal will run cmpress on startup."
                }
            }
        }
    }

    // 3. RNA-seq Active Alignment Source Verification
    if (!params.skip_rnaseq) {
        def has_samplesheet = params.rnaseq_samplesheet && file(params.rnaseq_samplesheet).exists()
        def has_bam = (params.bam && file(params.bam).exists()) || (params.bams)
        if (!has_samplesheet && !has_bam) {
            validation_errors.add("RNA-seq alignment is active, but neither a valid --rnaseq_samplesheet nor pre-computed --bam/--bams was provided.")
        }
    }

    // 4. Gene Prediction Reference Proteins non-empty check
    if (!params.skip_gene_prediction && !params.representative_gtf && params.proteins) {
        def prot_file = file(params.proteins)
        if (prot_file.isEmpty()) {
            validation_errors.add("Reference protein FASTA is empty: ${params.proteins}")
        }
    }

    // 5. Functional Annotation (InterProScan)
    if (!params.skip_functional_annotation && !params.interproscan_tsv && !params.interproscan_dir) {
        log.warn "[PRE-RUN GATE] --interproscan_dir not specified; functional annotation with InterProScan will be skipped."
    }

    // Halt immediately if any domain check failed
    if (validation_errors) {
        log.error "==================================================================="
        log.error " PRE-RUN VERIFICATION GATE FAILED: ${validation_errors.size()} CRITICAL CONFIGURATION ISSUE(S)"
        log.error "==================================================================="
        validation_errors.eachWithIndex { err, i ->
            log.error " [${i+1}] ${err}"
        }
        log.error "==================================================================="
        error "Pre-run verification gate failed. Please correct the parameters above before executing."
    }

    log.info "==================================================================="
    log.info " PRE-RUN VERIFICATION GATE PASSED: All inputs, databases, and paths verified."
    log.info "==================================================================="

    ch_assembly_file = file(params.assembly, checkIfExists: true)
    def clean_assembly_id = ch_assembly_file.name.replaceAll(/(\.gz|\.bgz)$/, '').replaceAll(/\.(fa|fasta|fna)$/, '')
    ch_assembly_meta = [ id: clean_assembly_id ]

    // Transparent on-the-fly decompression for (b)gzipped assemblies (.fa.gz, .fasta.gz, .fa.bgz)
    def is_compressed_assembly = params.assembly.toString().endsWith('.gz') || params.assembly.toString().endsWith('.bgz')
    if (is_compressed_assembly) {
        GUNZIP(
            channel.of([ ch_assembly_meta, ch_assembly_file ])
        )
        ch_assembly = GUNZIP.out.gunzip
    } else {
        ch_assembly = channel.of([ ch_assembly_meta, ch_assembly_file ])
    }

    // FASTA Index
    if (params.fai && !is_compressed_assembly) {
        ch_fai = channel.of(file(params.fai, checkIfExists: true))
    } else {
        SAMTOOLS_FAIDX(
            ch_assembly.map { meta, fasta -> [ meta, fasta, [] ] },
            false
        )
        ch_fai = SAMTOOLS_FAIDX.out.fai.map { _meta, fai -> fai }
    }

    // -----------------------------------------------------------------------
    // Stage 1: ncRNA Structural Annotation (Infernal cmscan)
    // Runs upfront to provide high-density tandem structural RNA masks (rDNA, 5S, snRNA arrays)
    // that shield Earl Grey / RepeatModeler from massive tandem unit explosions.
    // -----------------------------------------------------------------------
    ch_ncrna_gff = channel.empty()
    ch_structural_mask_bed = channel.empty()

    if (params.structural_mask_bed) {
        // Direct precomputed structural RNA tandem mask bypass
        ch_structural_mask_bed = channel.value([ ch_assembly_meta, file(params.structural_mask_bed, checkIfExists: true) ])
    }

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

    if (!params.structural_mask_bed && (!params.skip_ncrna || params.ncrna_gff)) {
        NCRNA_TANDEM_FILTER(
            ch_ncrna_gff.map { gff -> [ ch_assembly_meta, gff ] }
        )
        ch_structural_mask_bed = NCRNA_TANDEM_FILTER.out.bed
    }

    // -----------------------------------------------------------------------
    // Stage 2: Repeats & Satellites (Centroid Consensus Satellites + Earl Grey)
    // -----------------------------------------------------------------------
    ch_masked_assembly = ch_assembly
    ch_sat_masked_assembly = ch_assembly

    if (!params.skip_repeats) {
        REPEAT_ANNOTATION(
            ch_assembly,
            ch_fai,
            params.window_size,
            params.min_period,
            params.min_copies,
            params.min_identity,
            params.clade,
            params.skip_te,
            params.earlgrey_dir,
            params.earlgrey_seed_dir,
            ch_structural_mask_bed
        )
        ch_masked_assembly     = REPEAT_ANNOTATION.out.masked_fasta
        ch_sat_masked_assembly = REPEAT_ANNOTATION.out.sat_masked_fasta
    }

    // -----------------------------------------------------------------------
    // Stage 3: RNA-seq Alignment / Pre-computed Inspection
    // -----------------------------------------------------------------------
    ch_primary_bam = channel.empty()
    ch_individual_bams = channel.empty()
    ch_aligned_merged_bam = channel.empty()

    def precomputed_aln = params.bam ?: null

    if (precomputed_aln) {
        def aln_file = file(precomputed_aln, checkIfExists: true)
        def aln_idx = file("${aln_file}.bai").exists() ? file("${aln_file}.bai") :
                     (file("${aln_file.toString().replaceAll('\\.bam$', '')}.bai").exists() ? file("${aln_file.toString().replaceAll('\\.bam$', '')}.bai") :
                     (file("${aln_file}.crai").exists() ? file("${aln_file}.crai") : []))

        def ch_precomputed_raw = channel.of([ [ id: ch_assembly_meta.id ], aln_file, aln_idx ])
        def ch_aln_has_idx = ch_precomputed_raw.filter { _m, _bam, idx -> idx != [] }
        def ch_aln_needs_idx = ch_precomputed_raw.filter { _m, _bam, idx -> idx == [] }.map { m, bam, _idx -> [ m, bam ] }

        INDEX_PRECOMPUTED_BAM(ch_aln_needs_idx)
        def ch_aln_indexed = ch_aln_needs_idx.join(INDEX_PRECOMPUTED_BAM.out.index).map { m, bam, idx -> [ m, bam, idx ] }

        def ch_precomputed_with_idx = ch_aln_has_idx.mix(ch_aln_indexed)

        VALIDATE_RNASEQ_ALIGNMENT(
            ch_precomputed_with_idx,
            ch_fai.map { fai -> [ [ id: ch_assembly_meta.id ], fai ] }
        )

        // Evaluate whether alignment header has HISAT2 provenance and coordinate validity
        ch_checked_aln = VALIDATE_RNASEQ_ALIGNMENT.out.evaluation
        ch_hisat2_valid = ch_checked_aln.filter { _meta, _file, _idx, is_hisat2, is_coord ->
            is_hisat2.toString().trim() == 'true' && is_coord.toString().trim() == 'true'
        }

        ch_aligned_merged_bam = ch_hisat2_valid.map { meta, aln_f, _idx, _res, _coord -> [ meta, aln_f ] }
        ch_primary_bam = ch_hisat2_valid.map { _meta, aln_f, _idx, _res, _coord -> aln_f }
    }

    if (params.rnaseq_samplesheet) {
        ch_samplesheet = channel.fromPath(params.rnaseq_samplesheet, checkIfExists: true)

        // If no valid precomputed HISAT2 alignment was provided, run RNASEQ_ALIGN
        if (!precomputed_aln) {
            RNASEQ_ALIGN(
                ch_assembly,
                ch_fai,
                ch_samplesheet
            )
            ch_aligned_merged_bam = RNASEQ_ALIGN.out.merged_bam
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
        }
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
                : [ file(raw_bams_arg, checkIfExists: true) ]

            def ch_raw_sl_bams = channel.fromList(bam_file_list)
                .map { bfile ->
                    def meta = [ id: bfile.baseName.replaceAll('\\.slidr$', '').replaceAll('\\.sorted$', '') ]
                    def bai = file("${bfile}.bai").exists() ? file("${bfile}.bai") :
                             (file("${bfile.toString().replaceAll('\\.bam$', '')}.bai").exists() ? file("${bfile.toString().replaceAll('\\.bam$', '')}.bai") :
                             (file("${bfile}.crai").exists() ? file("${bfile}.crai") : []))
                    [ meta, bfile, bai ]
                }

            def ch_sl_has_idx = ch_raw_sl_bams.filter { _m, _bam, idx -> idx != [] }
            def ch_sl_needs_idx = ch_raw_sl_bams.filter { _m, _bam, idx -> idx == [] }.map { m, bam, _idx -> [ m, bam ] }

            INDEX_SL_BAM(ch_sl_needs_idx)
            def ch_sl_indexed = ch_sl_needs_idx.join(INDEX_SL_BAM.out.index).map { m, bam, idx -> [ m, bam, idx ] }

            ch_sl_bams = ch_sl_has_idx.mix(ch_sl_indexed)
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
    // Stage 5: Dual-Evidence Gene Prediction (or Direct Representative Bypass)
    // -----------------------------------------------------------------------
    ch_representative_gtf = channel.empty()
    ch_representative_faa = channel.empty()
    ch_representative_fna = channel.empty()

    if (params.representative_gtf) {
        file_rep_gtf = file(params.representative_gtf, checkIfExists: true)
        ch_representative_gtf = channel.of([ [ id: ch_assembly_meta.id ], file_rep_gtf ])

        // Expose/publish the representative GTF into standard output location
        ch_assembly.map { meta, _f ->
            def outdir_annot = file("${params.outdir}/annotation")
            outdir_annot.mkdirs()
            def target_gtf = file("${outdir_annot}/${meta.id}.representative.gtf")
            if (!target_gtf.exists()) {
                file_rep_gtf.copyTo(target_gtf)
            }
        }

        ch_plain_fasta = ch_assembly.map { _meta, fasta -> fasta }.first()

        if (params.representative_faa) {
            file_rep_faa = file(params.representative_faa, checkIfExists: true)
            ch_representative_faa = channel.of([ [ id: ch_assembly_meta.id ], file_rep_faa ])
        } else {
            GFFREAD_PROTEIN(
                ch_representative_gtf,
                ch_plain_fasta
            )
            ch_representative_faa = GFFREAD_PROTEIN.out.aa
        }

        GFFREAD_CDS(
            ch_representative_gtf,
            ch_plain_fasta
        )
        ch_representative_fna = GFFREAD_CDS.out.cds

        // -------------------------------------------------------------------
        // Stage 6: Functional Annotation (InterProScan)
        // -------------------------------------------------------------------
        if (!params.skip_functional_annotation) {
            if (params.interproscan_tsv) {
                file_ipr_tsv  = file(params.interproscan_tsv, checkIfExists: true)
                file_ipr_gff3 = params.interproscan_gff3 ? file(params.interproscan_gff3, checkIfExists: true) : null
                ch_assembly.map { meta, _f ->
                    def outdir_ipr = file("${params.outdir}/functional_annotation/interproscan")
                    outdir_ipr.mkdirs()
                    def target_ext_tsv = file_ipr_tsv.name.endsWith('.gz') ? '.iprscan.tsv.gz' : '.iprscan.tsv'
                    def target_tsv = file("${outdir_ipr}/${meta.id}${target_ext_tsv}")
                    if (!target_tsv.exists()) {
                        file_ipr_tsv.copyTo(target_tsv)
                    }
                    if (file_ipr_gff3) {
                        def target_ext_gff = file_ipr_gff3.name.endsWith('.gz') ? '.iprscan.gff3.gz' : '.iprscan.gff3'
                        def target_gff3 = file("${outdir_ipr}/${meta.id}${target_ext_gff}")
                        if (!target_gff3.exists()) {
                            file_ipr_gff3.copyTo(target_gff3)
                        }
                    }
                }
            } else if (params.interproscan_dir) {
                ch_ipr_dir = file(params.interproscan_dir, checkIfExists: true)
                FUNCTIONAL_ANNOTATION(
                    ch_representative_faa,
                    ch_ipr_dir
                )
            }
        }
    } else if (!params.skip_gene_prediction) {
        ch_proteins       = params.proteins ? file(params.proteins, checkIfExists: true) : []
        ch_busco_lineages = params.busco_lineages_path ? file(params.busco_lineages_path, checkIfExists: true) : []

        file_braker_gtf   = params.braker_gtf   ? file(params.braker_gtf, checkIfExists: true)   : null
        file_braker_aa    = params.braker_aa    ? file(params.braker_aa, checkIfExists: true)    : null
        file_braker_hints = params.braker_hints ? file(params.braker_hints, checkIfExists: true) : null
        file_helixer_gff3 = params.helixer_gff3 ? file(params.helixer_gff3, checkIfExists: true) : null
        file_helixer_gtf  = params.helixer_gtf  ? file(params.helixer_gtf, checkIfExists: true)  : null

        def skip_helixer_val = params.skip_helixer || !params.run_helixer

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
            file_helixer_gtf,
            skip_helixer_val
        )
        ch_representative_gtf = GENE_PREDICTION.out.representative_gtf
        ch_representative_faa = GENE_PREDICTION.out.representative_faa
        ch_representative_fna = GENE_PREDICTION.out.representative_fna

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
                    // Create convenient symlinks/copies in standard output location preserving compression
                    def target_ext_tsv = file_ipr_tsv.name.endsWith('.gz') ? '.iprscan.tsv.gz' : '.iprscan.tsv'
                    def target_tsv = file("${outdir_ipr}/${meta.id}${target_ext_tsv}")
                    if (!target_tsv.exists()) {
                        file_ipr_tsv.copyTo(target_tsv)
                    }
                    if (file_ipr_gff3) {
                        def target_ext_gff = file_ipr_gff3.name.endsWith('.gz') ? '.iprscan.gff3.gz' : '.iprscan.gff3'
                        def target_gff3 = file("${outdir_ipr}/${meta.id}${target_ext_gff}")
                        if (!target_gff3.exists()) {
                            file_ipr_gff3.copyTo(target_gff3)
                        }
                    }
                }
            } else if (params.interproscan_dir) {
                ch_ipr_dir = file(params.interproscan_dir, checkIfExists: true)
                ch_proteins_to_annotate = ch_representative_faa

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
    if (!params.skip_sl_discovery && (params.representative_gtf || !params.skip_gene_prediction)) {
        ch_first_rep_gtf = ch_representative_gtf.map { _meta, gtf -> gtf }.first()

        SLOPPR_OPERONS(
            SPLICED_LEADER_DISCOVERY.out.sl_bed,
            ch_first_rep_gtf,
            params.max_dist
        )
        ch_operons_gff3 = SLOPPR_OPERONS.out.operons_gff
    }

    // -----------------------------------------------------------------------
    // Stage 8: Transposable Element (TE) Multi-Evidence Scoring
    // -----------------------------------------------------------------------
    if ((params.representative_gtf || !params.skip_gene_prediction) && !params.skip_te_scoring) {
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
        ch_scorer_input = ch_representative_gtf
            .join(ch_representative_faa)
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
    if ((params.rnaseq_samplesheet || params.bam) && (params.representative_gtf || !params.skip_gene_prediction)) {
        ch_counts_input = ch_aligned_merged_bam
            .map { _meta, bam -> [ ch_assembly_meta.id, bam ] }
            .join(ch_representative_gtf.map { meta, gtf -> [ meta.id, gtf ] })
            .map { id, bam, gtf -> [ [ id: id, single_end: false ], [ bam ], gtf ] }

        SUBREAD_FEATURECOUNTS(
            ch_counts_input
        )
    }

    // -----------------------------------------------------------------------
    // Stage 10: Collate & Publish Pipeline Software Versions & Provenance
    // -----------------------------------------------------------------------
    def pipeline_version_str = workflow.commitId ? "${workflow.manifest.version} (commit ${workflow.commitId.take(9)})" : "${workflow.manifest.version}"
    def pipeline_prov_yaml = "Workflow:\n  nf-core/nematodeannotation: ${pipeline_version_str}"

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(channel.empty().mix(topic_versions.versions_file))
        .mix(channel.of(pipeline_prov_yaml))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true
        )
}
