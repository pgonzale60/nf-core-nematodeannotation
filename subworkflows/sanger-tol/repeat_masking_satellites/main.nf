//
// Subworkflow: REPEAT_MASKING_SATELLITES
// Implements De Novo Centroid Consensus Masking:
// 1. Chunked window harvesting with TRF v4.09 (bounded memory, zero OOM)
// 2. Filter repeats (period >= 50 bp, copies >= 100, identity >= 85%)
// 3. Cluster filtered units with VSEARCH at 90% identity to derive family centroids
// 4. Map centroids back to the assembly with BLASTN (E <= 0.01)
// 5. Merge adjacent intervals with bedtools merge -> non-overlapping satellite mask
// 6. Soft-mask assembly using bedtools maskfasta
//

include { BEDTOOLS_MAKEWINDOWS  } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BEDTOOLS_GETFASTA     } from '../../../modules/nf-core/bedtools/getfasta/main'
include { TRF                   } from '../../../modules/sanger-tol/trf/main'
include { SATELLITEFILTER       } from '../../../modules/sanger-tol/satellitefilter/main'
include { VSEARCH_CLUSTER       } from '../../../modules/nf-core/vsearch/cluster/main'
include { BLAST_MAKEBLASTDB     } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN          } from '../../../modules/sanger-tol/blast/blastn/main'
include { BLASTTOBED            } from '../../../modules/sanger-tol/blasttobed/main'
include { SATELLITEPROFILER     } from '../../../modules/sanger-tol/satelliteprofiler/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_SATELLITES } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MASKFASTA as BEDTOOLS_MASK_SATELLITES     } from '../../../modules/nf-core/bedtools/maskfasta/main'
include { BEDTOOLS_MASKFASTA as BEDTOOLS_HARDMASK_SATELLITES } from '../../../modules/nf-core/bedtools/maskfasta/main'

workflow REPEAT_MASKING_SATELLITES {
    take:
    ch_fasta           // channel: [ val(meta), path(fasta) ]
    ch_fai             // channel: path(fai) / genome file for windows
    _val_window_size   // int: window size for TRF harvesting (e.g. 500000)
    val_min_period     // int: min repeat period bp (default: 50)
    val_min_copies     // int/float: min copy number (default: 100)
    val_min_identity   // int: min sequence identity percent (default: 85)

    main:

    // 1. Generate non-overlapping genome windows to bound TRF memory
    ch_fai_meta = ch_fasta.map { meta, _fa -> meta }.first()
    ch_makewindows_input = ch_fai_meta.combine(ch_fai).map { meta, fai -> [ meta, fai ] }

    BEDTOOLS_MAKEWINDOWS(
        ch_makewindows_input
    )

    // 2. Split window BED into parallel chunks (e.g. 20 intervals per chunk)
    ch_bed_chunks = BEDTOOLS_MAKEWINDOWS.out.bed
        .flatMap { meta, bed ->
            bed.splitText(by: 20, file: true).collect { chunk_file ->
                [ [ id: "${meta.id}_chunk_${chunk_file.baseName}" ], chunk_file ]
            }
        }

    BEDTOOLS_GETFASTA(
        ch_bed_chunks,
        ch_fasta.map { _meta, fa -> fa }.first()
    )

    // 3. Parallel harvesting using TRF v4.09 across all chunks
    // Parameters: match=2, mismatch=7, delta=7, PM=80, PI=10, minscore=50, maxperiod=500
    TRF(
        BEDTOOLS_GETFASTA.out.fasta,
        2, 7, 7, 80, 10, 50, 500
    )

    // 4. Collect all chunk .dat files and filter repeats (period >= 50, copies >= 100, identity >= 85%)
    // collect() emits a flat List; combine() with a value channel re-flattens it into positional args.
    // Wrap in another list so combine sees a single tuple element: [meta, [dat1, dat2, ...]]
    ch_collected_dats = TRF.out.dat
        .map { _meta, dat -> dat }
        .collect()
        .map { dats -> [ dats ] }

    ch_filter_input = ch_fai_meta
        .combine(ch_collected_dats)
        .map { meta, dats -> [ meta, dats ] }

    SATELLITEFILTER(
        ch_filter_input,
        val_min_period,
        val_min_copies,
        val_min_identity
    )

    // 5. Cluster candidate consensus units with VSEARCH at 90% identity
    // Configure VSEARCH to use --cluster_fast and emit --centroids
    VSEARCH_CLUSTER(
        SATELLITEFILTER.out.fasta
    )

    // 6. Build BLAST database for the reference assembly
    BLAST_MAKEBLASTDB(
        ch_fasta,
        []
    )

    // 7. Map centroids back to the assembly with BLASTN (E <= 0.01)
    ch_centroids = VSEARCH_CLUSTER.out.centroids
        .map { meta, fasta -> [ [ id: "${meta.id}_centroids" ], fasta ] }

    ch_blast_db = BLAST_MAKEBLASTDB.out.db
        .map { meta, db -> [ meta, db ] }
        .first()

    BLAST_BLASTN(
        ch_centroids,
        ch_blast_db,
        [],
        [],
        []
    )

    // 8. Convert BLAST alignments to sorted BED intervals
    BLASTTOBED(
        BLAST_BLASTN.out.txtgz
    )

    // 9. Profile, classify, and annotate satellite arrays and monomers
    ch_profiler_input = BLASTTOBED.out.bed
        .join(ch_centroids)
        .map { meta, bed, centroids -> [ meta, bed, centroids ] }

    SATELLITEPROFILER(
        ch_profiler_input,
        100
    )

    // 10. Merge adjacent intervals into a non-overlapping satellite mask
    // Restore base assembly metadata so downstream processes retain assembly prefix
    ch_merged_sat_input = BLASTTOBED.out.bed
        .combine(ch_fasta.map { meta, _fa -> meta })
        .map { _bed_meta, bed, asm_meta -> [ asm_meta, bed ] }
    BEDTOOLS_MERGE_SATELLITES(
        ch_merged_sat_input
    )

    // 11. Mask assembly using bedtools maskfasta
    // 11a. Soft-masking (for downstream gene callers)
    BEDTOOLS_MASK_SATELLITES(
        BEDTOOLS_MERGE_SATELLITES.out.bed,
        ch_fasta.map { _meta, fa -> fa }.first()
    )

    // 11b. Hard-masking with Ns (strictly for Earl Grey / RepeatModeler de novo TE discovery)
    BEDTOOLS_HARDMASK_SATELLITES(
        BEDTOOLS_MERGE_SATELLITES.out.bed,
        ch_fasta.map { _meta, fa -> fa }.first()
    )

    emit:
    masked_fasta       = BEDTOOLS_MASK_SATELLITES.out.fasta       // [ meta, satellite_softmasked.fa ]
    hard_masked_fasta  = BEDTOOLS_HARDMASK_SATELLITES.out.fasta   // [ meta, satellite_hardmasked.fa ]
    satellite_bed      = BEDTOOLS_MERGE_SATELLITES.out.bed        // [ meta, aps4_satellites.bed ]
    centroids          = VSEARCH_CLUSTER.out.centroids            // [ meta, centroids.fasta.gz ]
    arrays_gff3        = SATELLITEPROFILER.out.arrays_gff3        // [ meta, *_satellite_arrays.gff3 ]
    monomers_gff3      = SATELLITEPROFILER.out.monomers_gff3      // [ meta, *_satellite_monomers.gff3 ]
    families_tsv       = SATELLITEPROFILER.out.families_tsv       // [ meta, *_satellite_families.tsv ]
    library_fasta      = SATELLITEPROFILER.out.library_fasta      // [ meta, *_satellite_monomers.fasta ]
    summary            = SATELLITEFILTER.out.summary              // [ meta, trf_summary.tsv ]
}
