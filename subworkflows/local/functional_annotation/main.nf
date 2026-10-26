//
// Subworkflow: FUNCTIONAL_ANNOTATION
// Automated functional annotation of predicted proteomes via InterProScan
//

include { INTERPROSCAN } from '../../../modules/local/interproscan/main'

workflow FUNCTIONAL_ANNOTATION {
    take:
    ch_proteins     // channel: [ val(meta), path(fasta) ]
    ch_ipr_dir      // path(interproscan_dir)

    main:

    INTERPROSCAN(
        ch_proteins,
        ch_ipr_dir
    )

    emit:
    tsv      = INTERPROSCAN.out.tsv
    gff3     = INTERPROSCAN.out.gff3
    xml      = INTERPROSCAN.out.xml
}
