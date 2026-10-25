//
// Subworkflow: REPEAT_ANNOTATION
// High-performance repeat identification, satellite profiling, and transposable element masking.
//
// Protocol sequence:
// 1. REPEAT_MASKING_SATELLITES:
//    - Chunked window TRF v4.09 harvesting (memory bounded)
//    - High-copy filter (period >= 50 bp, copies >= 100, identity >= 85%)
//    - VSEARCH centroid clustering (90% identity)
//    - Assembly-wide BLASTN mapping
//    - SATELLITEPROFILER: canonical circular phasing, monomer parsing, array metrics
//    - Soft-masked and hard-masked satellite FASTA generation
// 2. EARLGREY:
//    - Interspersed TE discovery on satellite hard-masked assembly (prevents TE callers
//      from misclassifying multimegabase satellite arrays as novel TEs)
//    - Produces final fully soft-masked FASTA and TE GFF3/BED/library
//

include { REPEAT_MASKING_SATELLITES } from '../repeat_masking_satellites/main'
include { EARLGREY                  } from '../../../modules/sanger-tol/earlgrey/main'
include { BEDTOOLS_MASKFASTA as BEDTOOLS_MASK_TE } from '../../../modules/nf-core/bedtools/maskfasta/main'
include { BEDTOOLS_MASKFASTA as BEDTOOLS_MASK_STRUCTURAL } from '../../../modules/nf-core/bedtools/maskfasta/main'

workflow REPEAT_ANNOTATION {
    take:
    ch_fasta               // channel: [ val(meta), path(fasta) ]
    ch_fai                 // channel: path(fai)
    val_window_size        // int: window size for TRF (default: 500000)
    val_min_period         // int: min period (default: 50)
    val_min_copies         // int/float: min copies (default: 100)
    val_min_identity       // int: min identity (default: 85)
    val_rhabditida_clade   // string: Earl Grey clade (default: 'rhabditida')
    val_skip_te            // boolean: skip Earl Grey (default: false)
    val_earlgrey_dir       // string/path: optional precomputed Earl Grey directory (default: null)
    ch_structural_mask_bed // channel: optional [ val(meta), path(bed) ] of structural RNA tandem arrays to hard-mask prior to Earl Grey

    main:

    // 1. Satellite Consensus Masking & Profiling (Stevens protocol + SatelliteProfiler)
    REPEAT_MASKING_SATELLITES(
        ch_fasta,
        ch_fai,
        val_window_size,
        val_min_period,
        val_min_copies,
        val_min_identity
    )

    ch_sat_masked_fasta      = REPEAT_MASKING_SATELLITES.out.masked_fasta
    ch_sat_hard_masked_fasta = REPEAT_MASKING_SATELLITES.out.hard_masked_fasta

    // 1b. Hard-mask structural RNA tandem arrays (rDNA, 5S, snRNA arrays) onto sat_hard_masked_fasta
    // to shield Earl Grey / RepeatModeler from massive tandem consensus explosions.
    ch_earlgrey_input_fasta = ch_sat_hard_masked_fasta
    if (ch_structural_mask_bed) {
        BEDTOOLS_MASK_STRUCTURAL(
            ch_structural_mask_bed,
            ch_sat_hard_masked_fasta.map { _meta, fa -> fa }
        )
        ch_earlgrey_input_fasta = BEDTOOLS_MASK_STRUCTURAL.out.fasta
    }
    ch_sat_bed               = REPEAT_MASKING_SATELLITES.out.satellite_bed
    ch_centroids             = REPEAT_MASKING_SATELLITES.out.centroids
    ch_arrays_gff3           = REPEAT_MASKING_SATELLITES.out.arrays_gff3
    ch_monomers_gff3         = REPEAT_MASKING_SATELLITES.out.monomers_gff3
    ch_families_tsv          = REPEAT_MASKING_SATELLITES.out.families_tsv
    ch_library_fasta         = REPEAT_MASKING_SATELLITES.out.library_fasta

    // 2. Interspersed Transposable Element Masking (Earl Grey or precomputed deliverables)
    ch_fully_masked_fasta = ch_sat_masked_fasta
    ch_te_gff3            = channel.empty()
    ch_te_bed             = channel.empty()

    if (val_earlgrey_dir) {
        // Ingest precomputed Earl Grey deliverables directly
        ch_precomputed_eg = ch_fasta.map { item ->
            def meta = item[0]
            def eg_dir = file(val_earlgrey_dir, checkIfExists: true)
            def bed_files = file("${eg_dir}/*_summaryFiles/*.filteredRepeats.bed")
            def gff_files = file("${eg_dir}/*_summaryFiles/*.filteredRepeats.gff*")
            def te_bed = bed_files ? bed_files[0] : file("${eg_dir}/*.filteredRepeats.bed", checkIfExists: true)
            def te_gff = gff_files ? gff_files[0] : file("${eg_dir}/*.filteredRepeats.gff*", checkIfExists: true)
            return [ meta, te_bed, te_gff ]
        }

        ch_te_bed  = ch_precomputed_eg.map { item -> [ item[0], item[1] ] }
        ch_te_gff3 = ch_precomputed_eg.map { item -> [ item[0], item[2] ] }

        BEDTOOLS_MASK_TE(
            ch_te_bed,
            ch_sat_masked_fasta.map { _meta, fa -> fa }
        )
        ch_fully_masked_fasta = BEDTOOLS_MASK_TE.out.fasta

    } else if (!val_skip_te) {
        EARLGREY(
            ch_earlgrey_input_fasta,
            val_rhabditida_clade
        )
        ch_te_gff3  = EARLGREY.out.gff3
        ch_te_bed   = EARLGREY.out.bed

        // Soft-mask the Earl Grey TE BED directly onto the satellite soft-masked FASTA.
        // This ensures the final assembly is completely soft-masked (both satellites and TEs in lowercase),
        // with ZERO hardmasked Ns introduced into the assembly.
        BEDTOOLS_MASK_TE(
            ch_te_bed,
            ch_sat_masked_fasta.map { _meta, fa -> fa }
        )
        ch_fully_masked_fasta = BEDTOOLS_MASK_TE.out.fasta
    }

    emit:
    masked_fasta          = ch_fully_masked_fasta          // [ meta, fully_softmasked.fa ]
    sat_masked_fasta      = ch_sat_masked_fasta            // [ meta, sat_softmasked.fa ]
    hard_masked_fasta     = ch_sat_hard_masked_fasta       // [ meta, sat_hardmasked.fa ]
    satellite_bed         = ch_sat_bed                     // [ meta, satellites.bed ]
    satellite_centroids   = ch_centroids                   // [ meta, centroids.fasta.gz ]
    satellite_arrays_gff3 = ch_arrays_gff3                 // [ meta, *_satellite_arrays.gff3 ]
    satellite_monomers_gff3 = ch_monomers_gff3             // [ meta, *_satellite_monomers.gff3 ]
    satellite_families_tsv = ch_families_tsv              // [ meta, *_satellite_families.tsv ]
    satellite_library_fasta = ch_library_fasta             // [ meta, *_satellite_monomers.fasta ]
    te_gff3               = ch_te_gff3                     // [ meta, te.gff3 ]
    te_bed                = ch_te_bed                      // [ meta, te.bed ]
}
