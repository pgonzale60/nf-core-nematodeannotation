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
include { CHECK_HISAT2_ALIGNMENT    } from './modules/local/rnaseq_check/main'
include { GFFREAD_PROTEIN           } from './modules/sanger-tol/gffread/main'
include { GFFREAD_CDS               } from './modules/sanger-tol/gffread/main'
include { softwareVersionsToYAML    } from './subworkflows/nf-core/utils_nfcore_pipeline/main'

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
        log.info """
        ===================================================================
        nf-core-nematodeannotation: Nematode Genome Annotation Pipeline
        ===================================================================
        Usage:
          nextflow run main.nf --assembly <path/to/assembly.fa> [options]

        Mandatory arguments:
          --assembly             Path to input FASTA file (.fa, .fasta, .fa.gz, .fasta.gz)

        Input / Output options:
          --outdir               Path to directory where results are published
          --fai                  Path to assembly FASTA index (.fai) [optional]
          --rnaseq_samplesheet   Path to RNA-seq samplesheet CSV [optional]
          --proteins             Path to OrthoDB/lineage protein FASTA for BRAKER3 [optional]
          --busco_lineage        BUSCO lineage dataset (default: 'chromadorea_odb12')
          --interproscan_dir     Path to InterProScan data directory [optional]

        Bypass / Pre-computed options:
          --representative_gtf   Pre-computed winning gene models GTF (bypasses BRAKER/Helixer/TSEBRA tournament)
          --representative_faa   Pre-computed winning proteins FASTA [optional, derived from GTF if omitted]
          --braker_gtf           Pre-computed BRAKER GTF
          --braker_aa            Pre-computed BRAKER proteins
          --helixer_gff3         Pre-computed Helixer GFF3
          --interproscan_tsv     Pre-computed InterProScan TSV
          --satellites_bed       Pre-computed satellites BED
          --rfam_cm              Rfam covariance models file
        ===================================================================
        """.stripIndent()
        exit 0
    }

    def commit_suffix = workflow.commitId ? workflow.commitId.take(9) : 'dev'

    if (!params.assembly) {
        error "Mandatory parameter --assembly <path/to/assembly.fa> not specified. Run with --help for usage."
    }

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
    // Stage 1: Repeats & Satellites (Centroid Consensus Satellites + Earl Grey)
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
    // Stage 3: RNA-seq Alignment / Pre-computed Inspection
    // -----------------------------------------------------------------------
    ch_primary_bam = channel.empty()
    ch_individual_bams = channel.empty()
    ch_aligned_merged_bam = channel.empty()

    def precomputed_aln = params.bam ?: null
    def needs_alignment = true

    if (precomputed_aln) {
        def aln_file = file(precomputed_aln, checkIfExists: true)
        def aln_idx = file("${aln_file}.bai").exists() ? file("${aln_file}.bai") :
                     (file("${aln_file}.crai").exists() ? file("${aln_file}.crai") : [])

        CHECK_HISAT2_ALIGNMENT(
            channel.of([ [ id: ch_assembly_meta.id ], aln_file, aln_idx ]),
            ch_assembly.map { _meta, fasta -> fasta }.first()
        )

        // Evaluate whether header has HISAT2 provenance
        // In DSL2 Nextflow scripts, we can check or branch
        ch_checked_aln = CHECK_HISAT2_ALIGNMENT.out.evaluation
        ch_hisat2_valid = ch_checked_aln.filter { _meta, _file, is_hisat2 -> is_hisat2.toString().trim() == 'true' }
        ch_hisat2_invalid = ch_checked_aln.filter { _meta, _file, is_hisat2 -> is_hisat2.toString().trim() != 'true' }

        ch_aligned_merged_bam = ch_hisat2_valid.map { meta, aln_f, _res -> [ meta, aln_f ] }
        ch_primary_bam = ch_hisat2_valid.map { _meta, aln_f, _res -> aln_f }
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
            .map { meta, bam -> [ ch_assembly_meta.id, bam ] }
            .join(ch_representative_gtf.map { meta, gtf -> [ meta.id, gtf ] })
            .map { id, bam, gtf -> [ [ id: id, single_end: false ], [ bam ], gtf ] }

        SUBREAD_FEATURECOUNTS(
            ch_counts_input
        )
    }

    // -----------------------------------------------------------------------
    // Stage 10: Collate & Publish Pipeline Software Versions
    // -----------------------------------------------------------------------
    def topic_versions = Channel.topic("versions")
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
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true
        )
}
