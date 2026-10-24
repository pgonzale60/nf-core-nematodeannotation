process TESCORER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/python:3.10.4"

    input:
    tuple val(meta), path(gtf), path(fasta), path(iprscan_tsv), path(earlgrey_gff3), path(satellites_bed), path(operons_gff3), path(orthogroups_tsv)

    output:
    tuple val(meta), path("${prefix}.gene_te_scores.tsv"), emit: te_scores_tsv
    tuple val(meta), path("${prefix}.te_scoring_summary.md"), emit: summary_md
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def ipr_arg    = (iprscan_tsv && !(iprscan_tsv.name == 'NO_FILE')) ? "--iprscan ${iprscan_tsv}" : ""
    def eg_arg     = (earlgrey_gff3 && !(earlgrey_gff3.name == 'NO_FILE')) ? "--earlgrey ${earlgrey_gff3}" : ""
    def sat_arg    = (satellites_bed && !(satellites_bed.name == 'NO_FILE')) ? "--satellites ${satellites_bed}" : ""
    def op_arg     = (operons_gff3 && !(operons_gff3.name == 'NO_FILE')) ? "--operons ${operons_gff3}" : ""
    def og_arg     = (orthogroups_tsv && !(orthogroups_tsv.name == 'NO_FILE')) ? "--orthogroups ${orthogroups_tsv}" : ""
    """
    # Tool: prep_longest_orthofinder.py | Script-SHA: ee56dbe68ae840b7
    python3 ${moduleDir}/resources/usr/bin/te_scorer.py \\
        --gtf "${gtf}" \\
        --fasta "${fasta}" \\
        ${ipr_arg} \\
        ${eg_arg} \\
        ${sat_arg} \\
        ${op_arg} \\
        ${og_arg} \\
        --out-tsv "${prefix}.gene_te_scores.tsv" \\
        --summary "${prefix}.te_scoring_summary.md"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.gene_te_scores.tsv"
    touch "${prefix}.te_scoring_summary.md"
    """
}
