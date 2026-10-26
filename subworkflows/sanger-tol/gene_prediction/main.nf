//
// Subworkflow: GENE_PREDICTION
// Multi-Predictor Tournament (BRAKER3, Augustus.hints, GeneMark, Helixer, TSEBRA consensus)
// Evaluated via consolidated in-process protein BUSCO to select the optimal representative model.
//

include { BRAKER3                                       } from '../../../modules/sanger-tol/braker3/main'
include { HELIXER                                       } from '../../../modules/sanger-tol/helixer/main'
include { GFFREAD_HELIXER                               } from '../../../modules/sanger-tol/gffread/main'
include { GFFREAD_PROTEIN as GFFREAD_GENEMARK_PROTEIN   } from '../../../modules/sanger-tol/gffread/main'
include { GFFREAD_PROTEIN as GFFREAD_HELIXER_PROTEIN    } from '../../../modules/sanger-tol/gffread/main'
include { GFFREAD_PROTEIN as GFFREAD_TSEBRA_PROTEIN     } from '../../../modules/sanger-tol/gffread/main'
include { GFFREAD_PROTEIN as GFFREAD_REPRESENTATIVE_FAA } from '../../../modules/sanger-tol/gffread/main'
include { GFFREAD_CDS as GFFREAD_REPRESENTATIVE_FNA     } from '../../../modules/sanger-tol/gffread/main'
include { TSEBRA                                        } from '../../../modules/sanger-tol/tsebra/main'
include { SELECTBESTMODEL                               } from '../../../modules/sanger-tol/selectbestmodel/main'
include { ARCHIVE_SUBOPTIMAL_MODELS                     } from '../../../modules/local/suboptimal_models/main'

workflow GENE_PREDICTION {
    take:
    ch_assembly            // channel: [ val(meta), path(fasta) ]
    ch_bam                 // channel: path(bam) or []
    ch_proteins            // channel: path(fasta) or []
    ch_busco_lineages      // channel: path(lineages_dir) or []
    val_helixer_lineage    // string: Helixer lineage model (e.g. 'invertebrate')
    val_busco_lineage      // string: BUSCO lineage dataset (e.g. 'chromadorea_odb12')
    file_braker_gtf        // path or null: pre-computed BRAKER3 GTF
    file_braker_aa         // path or null: pre-computed BRAKER3 proteins
    file_braker_hints      // path or null: pre-computed hintsfile.gff
    file_helixer_gff3      // path or null: pre-computed Helixer GFF3
    file_helixer_gtf       // path or null: pre-computed Helixer GTF
    val_skip_helixer       // boolean: whether to skip Helixer execution (default: true)

    main:
    ch_assembly_plain = ch_assembly.map { _meta, fasta -> fasta }

    // 1. BRAKER3 (RNA-seq BAM + Protein Homology or pre-computed)
    if (file_braker_gtf && file_braker_aa) {
        ch_braker_gtf   = ch_assembly.map { meta, _f -> [ meta, file_braker_gtf ] }
        ch_braker_aa    = ch_assembly.map { meta, _f -> [ meta, file_braker_aa ] }
        ch_augustus_gtf = ch_assembly.map { meta, _f -> [ meta, [] ] }
        ch_augustus_aa  = ch_assembly.map { meta, _f -> [ meta, [] ] }
        ch_genemark_gtf = ch_assembly.map { meta, _f -> [ meta, [] ] }
        ch_braker_hints = file_braker_hints ? file_braker_hints : []
    } else {
        BRAKER3(
            ch_assembly,
            ch_bam,
            [],
            [],
            ch_proteins,
            []
        )
        ch_braker_gtf   = BRAKER3.out.gtf
        ch_braker_aa    = BRAKER3.out.aa
        ch_augustus_gtf = BRAKER3.out.augustus_gtf
        ch_augustus_aa  = BRAKER3.out.augustus_aa
        ch_genemark_gtf = BRAKER3.out.genemark_gtf
        ch_braker_hints = BRAKER3.out.hintsfile.map { _meta, hints -> hints }
    }

    // Translate GeneMark GTF to AA if available
    ch_genemark_for_prot = ch_genemark_gtf.map { meta, gtf -> [ [ id: "${meta.id}.genemark" ], gtf ] }
    GFFREAD_GENEMARK_PROTEIN(
        ch_genemark_for_prot,
        ch_assembly_plain
    )
    ch_genemark_aa = GFFREAD_GENEMARK_PROTEIN.out.aa.map { meta, aa -> [ [ id: meta.id.replaceAll(/\.genemark$/, '') ], aa ] }

    // 2. HELIXER (Deep Learning ab initio or pre-computed, optional)
    ch_helixer_gtf = channel.empty()
    ch_helixer_aa  = channel.empty()
    ch_cand_helixer = channel.empty()

    if (file_helixer_gtf) {
        ch_helixer_gtf = ch_assembly.map { meta, _f -> [ meta, file_helixer_gtf ] }
    } else if (file_helixer_gff3) {
        GFFREAD_HELIXER(
            ch_assembly.map { meta, _f -> [ meta, file_helixer_gff3 ] }
        )
        ch_helixer_gtf = GFFREAD_HELIXER.out.gtf
    } else if (!val_skip_helixer) {
        HELIXER(
            ch_assembly,
            val_helixer_lineage
        )

        // 3. Convert Helixer GFF3 to clean GTF
        GFFREAD_HELIXER(
            HELIXER.out.gff3
        )
        ch_helixer_gtf = GFFREAD_HELIXER.out.gtf
    }

    // 4. TSEBRA Consensus Filtering & Reconciliation
    // Combine BRAKER3 GTF with Helixer GTF if available, otherwise filter BRAKER3 models alone
    ch_combined_gtfs = ch_braker_gtf
        .join(ch_helixer_gtf, remainder: true)
        .map { meta, gtf1, gtf2 ->
            def gtfs = [ gtf1 ]
            if (gtf2) {
                gtfs.add(gtf2)
            }
            [ meta, gtfs ]
        }

    TSEBRA(
        ch_combined_gtfs,
        ch_braker_hints,
        [],
        []
    )

    // 5. Extract protein sequences for Helixer (if active) and TSEBRA consensus
    if (file_helixer_gtf || file_helixer_gff3 || !val_skip_helixer) {
        ch_helixer_gtf_for_prot = ch_helixer_gtf.map { meta, gtf -> [ [ id: "${meta.id}.helixer" ], gtf ] }
        GFFREAD_HELIXER_PROTEIN(
            ch_helixer_gtf_for_prot,
            ch_assembly_plain
        )
        ch_helixer_aa  = GFFREAD_HELIXER_PROTEIN.out.aa.map { meta, aa -> [ [ id: meta.id.replaceAll(/\.helixer$/, '') ], aa ] }

        ch_cand_helixer = ch_helixer_gtf.map { meta, gtf -> [ meta.id, meta, gtf ] }
            .join(ch_helixer_aa.map { meta, aa -> [ meta.id, aa ] })
            .map { id, meta, gtf, aa -> [ id, meta, 'Helixer', gtf, aa ] }
    }

    ch_tsebra_gtf_for_prot = TSEBRA.out.tsebra_gtf.map { meta, gtf -> [ [ id: "${meta.id}.tsebra" ], gtf ] }
    GFFREAD_TSEBRA_PROTEIN(
        ch_tsebra_gtf_for_prot,
        ch_assembly_plain
    )
    ch_tsebra_aa  = GFFREAD_TSEBRA_PROTEIN.out.aa.map { meta, aa -> [ [ id: meta.id.replaceAll(/\.tsebra$/, '') ], aa ] }

    // 6. Consolidated Tournament Selection via in-process BUSCO
    ch_cand_braker = ch_braker_gtf.map { meta, gtf -> [ meta.id, meta, gtf ] }
        .join(ch_braker_aa.map { meta, aa -> [ meta.id, aa ] })
        .map { id, meta, gtf, aa -> [ id, meta, 'BRAKER3', gtf, aa ] }

    ch_cand_augustus = ch_augustus_gtf.map { meta, gtf -> [ meta.id, meta, gtf ] }
        .join(ch_augustus_aa.map { meta, aa -> [ meta.id, aa ] })
        .map { id, meta, gtf, aa -> [ id, meta, 'Augustus_hints', gtf, aa ] }

    ch_cand_genemark = ch_genemark_gtf.map { meta, gtf -> [ meta.id, meta, gtf ] }
        .join(ch_genemark_aa.map { meta, aa -> [ meta.id, aa ] })
        .map { id, meta, gtf, aa -> [ id, meta, 'GeneMark', gtf, aa ] }

    ch_cand_tsebra = TSEBRA.out.tsebra_gtf.map { meta, gtf -> [ meta.id, meta, gtf ] }
        .join(ch_tsebra_aa.map { meta, aa -> [ meta.id, aa ] })
        .map { id, meta, gtf, aa -> [ id, meta, 'TSEBRA_Consensus', gtf, aa ] }

    ch_tournament_input = ch_cand_braker
        .mix(ch_cand_augustus, ch_cand_genemark, ch_cand_helixer, ch_cand_tsebra)
        .groupTuple(by: 0)
        .map { id, metas, labels, gtfs, aas ->
            [ metas[0], labels, gtfs, aas ]
        }

    SELECTBESTMODEL(
        ch_tournament_input,
        val_busco_lineage,
        ch_busco_lineages
    )

    // Standardize representative outputs to NCBI/EBI .faa and .fna
    GFFREAD_REPRESENTATIVE_FAA(
        SELECTBESTMODEL.out.representative_gtf,
        ch_assembly_plain
    )

    GFFREAD_REPRESENTATIVE_FNA(
        SELECTBESTMODEL.out.representative_gtf,
        ch_assembly_plain
    )

    // Archive all suboptimal GTF models into <id>.suboptimal_gene_models.tar.gz
    ch_all_gtfs = ch_braker_gtf.map { meta, gtf -> [ meta.id, meta, gtf ] }
        .mix(
            ch_augustus_gtf.map { meta, gtf -> [ meta.id, meta, gtf ] },
            ch_genemark_gtf.map { meta, gtf -> [ meta.id, meta, gtf ] },
            ch_helixer_gtf.map  { meta, gtf -> [ meta.id, meta, gtf ] },
            TSEBRA.out.tsebra_gtf.map { meta, gtf -> [ meta.id, meta, gtf ] }
        )
        .groupTuple(by: 0)
        .map { id, metas, gtfs -> [ metas[0], gtfs.findAll { it != [] } ] }

    ARCHIVE_SUBOPTIMAL_MODELS(
        ch_all_gtfs
    )

    emit:
    representative_gtf        = SELECTBESTMODEL.out.representative_gtf
    representative_faa        = GFFREAD_REPRESENTATIVE_FAA.out.aa
    representative_fna        = GFFREAD_REPRESENTATIVE_FNA.out.cds
    representative_tsv        = SELECTBESTMODEL.out.decision_tsv
    suboptimal_models_archive = ARCHIVE_SUBOPTIMAL_MODELS.out.archive
}
