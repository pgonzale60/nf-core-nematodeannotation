//
// Subworkflow: SPLICED_LEADER_DISCOVERY
// Upstream annotation-agnostic spliced leader discovery and TSS hint generation.
//
// Designed to run upstream of BRAKER or any other gene annotation pipeline:
//   1. SLIDR_EXTRACTTAILS: Extract 5' unaligned soft clips & trans-splice junctions per BAM
//   2. SLIDR_CLUSTERTAILS: Pool tails across libraries and cluster 3'-anchored
//   3. SLIDR_SLRNASITES:   Align centroids to assembly to detect potential leader gene loci
//   4. SLIDR_SLRNAFILTER:  Validate canonical donor (GT) and Sm-binding motifs de novo
//   5. SLIDR_CONSENSUS:    Attribute acceptor sites to discovered leaders and tally usage
//   6. SLIDR_TSSEXPORTER:  Export empirical TSS hints (GFF3) for downstream gene callers (BRAKER)
//
// NO gene annotation is required or assumed. Operon demarcation is deliberately omitted
// here so that TSS hints can directly inform primary gene predictor models.
//

include { SLIDR_EXTRACTTAILS } from '../../../modules/local/slidr/extracttails/main'
include { SLIDR_CLUSTERTAILS } from '../../../modules/local/slidr/clustertails/main'
include { SLIDR_SLRNASITES   } from '../../../modules/local/slidr/slrnasites/main'
include { SLIDR_SLRNAFILTER  } from '../../../modules/local/slidr/slrnafilter/main'
include { SLIDR_CONSENSUS    } from '../../../modules/local/slidr/consensus/main'
include { SLIDR_TSSEXPORTER  } from '../../../modules/local/slidr/tssexporter/main'

workflow SPLICED_LEADER_DISCOVERY {
    take:
    ch_bams                          // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta                         // channel: [ val(meta), path(fasta) ]
    ch_ncrna_gff                     // channel: path(gff3) from NCRNA_ANNOTATION or []
    val_min_support                  // int: min supporting reads per acceptor site (default: 5)
    val_min_tail_len                 // int: min soft-clipped tail length retained (default: 10)
    val_acceptor                     // str: acceptor site regex (default: AG)
    val_max_leader_gene_distance     // int: discard sites this close to leader gene (default: 1000)
    val_sl_donor                     // str: leader RNA splice donor regex (default: GT)
    val_sl_sm                        // str: Sm-binding-site regex (default: AATTTTTGGA|AATTTTGGA|AATTTGGA)
    val_sl_max_rna_len               // int: max leader RNA length excluding leader (default: 120)
    val_sl_max_overlap               // int: max leader overlap with donor (default: 3)
    val_sl_evalue                    // num: BLASTN e-value for centroid mapping (default: 10)
    val_min_tss_reads                // int: min pooled reads per site for TSS export (default: 3)
    val_min_tss_libraries            // int: min independent libraries per site (default: 1)
    val_rrna_flank                   // int: flank distance (bp) to exclude around rRNA arrays (default: 500)

    main:
    ch_rrna_gff = ch_ncrna_gff ? ch_ncrna_gff.first() : []

    // 1. Parallel soft-clipped tail and trans-splice acceptor extraction per BAM
    ch_extract_input = ch_bams
        .combine(ch_fasta)
        .map { meta_bam, bam, bai, _meta_fa, fasta ->
            [ meta_bam, bam, bai, fasta ]
        }

    SLIDR_EXTRACTTAILS(
        ch_extract_input.map { meta, bam, bai, _fa -> [ meta, bam, bai ] },
        ch_fasta.map { _meta, fa -> fa }.first()
    )

    // 2. Pool tails across libraries and cluster them with 3' anchoring
    ch_pooled_tails = SLIDR_EXTRACTTAILS.out.tails
        .map { _meta, tails -> tails }
        .collect()

    SLIDR_CLUSTERTAILS(
        ch_pooled_tails.map { files -> [ [ id: 'pooled' ], files ] },
        val_min_tail_len
    )

    // 3. Map cluster centroids to the genome to detect leader RNA genes
    SLIDR_SLRNASITES(
        SLIDR_CLUSTERTAILS.out.centroids,
        ch_fasta.first(),
        val_min_support,
        val_sl_evalue
    )

    // 4. Filter leader candidates by splice donor and Sm-binding site motif
    SLIDR_SLRNAFILTER(
        SLIDR_SLRNASITES.out.blast,
        SLIDR_CLUSTERTAILS.out.centroids,
        ch_fasta.first(),
        val_sl_donor,
        val_sl_sm,
        val_sl_max_rna_len,
        val_sl_max_overlap
    )

    // 5. Cross-reference acceptor sites against validated leaders, excluding rRNA mapping reads
    ch_pooled_sites = SLIDR_EXTRACTTAILS.out.sites
        .map { _meta, sites -> sites }
        .collect()
        .map { files -> [ [ id: 'pooled' ], files ] }
        .first()

    SLIDR_CONSENSUS(
        ch_pooled_sites,
        SLIDR_CLUSTERTAILS.out.readinfo.map { _meta, f -> f }.first(),
        SLIDR_SLRNAFILTER.out.tsv.map { _meta, f -> f }.first(),
        ch_rrna_gff,
        [],
        val_min_support,
        val_acceptor,
        val_max_leader_gene_distance,
        val_rrna_flank
    )

    // 6. Export high-confidence empirical TSS hints for gene annotation
    SLIDR_TSSEXPORTER(
        SLIDR_CONSENSUS.out.site_pool,
        val_min_tss_reads,
        val_min_tss_libraries
    )

    emit:
    tails_fa         = SLIDR_EXTRACTTAILS.out.tails        // [ meta, tails.fa.gz ]
    tail_stats       = SLIDR_EXTRACTTAILS.out.stats        // [ meta, tail_stats.json ]
    leader_centroids = SLIDR_CLUSTERTAILS.out.centroids     // [ meta, centroids.fa ]
    leader_blast     = SLIDR_SLRNASITES.out.blast           // [ meta, centroids.blastn.tsv.gz ]
    leader_rnas      = SLIDR_SLRNAFILTER.out.tsv           // [ meta, sl_rna_candidates.tsv ]
    leader_rna_gff   = SLIDR_SLRNAFILTER.out.gff         // [ meta, sl_rna_candidates.gff3 ]
    sl_bed           = SLIDR_CONSENSUS.out.bed             // [ meta, sl_acceptor_sites.bed ]
    sl_gff           = SLIDR_CONSENSUS.out.gff             // [ meta, sl_acceptor_sites.gff3 ]
    sl_summary       = SLIDR_CONSENSUS.out.summary         // [ meta, sl_summary.tsv ]
    sl_fasta         = SLIDR_CONSENSUS.out.sl_fasta        // [ meta, sl_sequences.fa ]
    tss_gff          = SLIDR_TSSEXPORTER.out.tss_gff       // [ meta, prefix.tss.gff3 ]
    tss_metrics      = SLIDR_TSSEXPORTER.out.metrics       // [ meta, prefix.tss_metrics.tsv ]
}
