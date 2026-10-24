//
// Subworkflow: NCRNA_ANNOTATION
// Non-coding RNA structural annotation with Infernal cmscan against Rfam
// and generation of an exclusion mask for downstream gene callers.
//

include { INFERNAL_CMSCAN                        } from '../../../modules/sanger-tol/infernal/cmscan/main'
include { INFERNAL_CMSCANTOGFF3                  } from '../../../modules/sanger-tol/infernal/cmscantogff3/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_NCRNA } from '../../../modules/nf-core/bedtools/merge/main'

workflow NCRNA_ANNOTATION {
    take:
    ch_fasta    // channel: [ val(meta), path(fasta) ]
    ch_rfam_cm  // channel: path(cmfile)

    main:
    ch_versions = channel.empty()

    // 1. Monolithic heavy compute search (pure engine)
    INFERNAL_CMSCAN(
        ch_fasta,
        ch_rfam_cm
    )
    ch_versions = ch_versions.mix(INFERNAL_CMSCAN.out.versions_infernal)

    // 2. Decoupled lightweight format conversion to GFF3 and raw exclusion BED
    INFERNAL_CMSCANTOGFF3(
        INFERNAL_CMSCAN.out.tbl
    )
    ch_versions = ch_versions.mix(INFERNAL_CMSCANTOGFF3.out.versions_perl)

    // 3. Merge overlapping / adjacent ncRNA features to construct gene caller exclusion BED
    BEDTOOLS_MERGE_NCRNA(
        INFERNAL_CMSCANTOGFF3.out.exclusion_bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE_NCRNA.out.versions_bedtools)

    emit:
    gff3          = INFERNAL_CMSCANTOGFF3.out.gff3    // [ meta, ncrna.gff3 ]
    tbl           = INFERNAL_CMSCAN.out.tbl           // [ meta, ncrna.tbl.gz ]
    exclusion_bed = BEDTOOLS_MERGE_NCRNA.out.bed      // [ meta, ncrna_exclusion.bed ]
    versions      = ch_versions
}
