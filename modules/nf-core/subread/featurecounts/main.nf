process SUBREAD_FEATURECOUNTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/subread:2.1.1--h577a1d6_0'
        : 'quay.io/biocontainers/subread:2.1.1--h577a1d6_0'}"

    input:
    tuple val(meta), path(bams), path(annotation)

    output:
    tuple val(meta), path("*featureCounts.tsv"), emit: counts
    tuple val(meta), path("*featureCounts.tsv.summary"), emit: summary
    tuple val("${task.process}"), val('subread'), eval("featureCounts -v 2>&1 | sed 's/featureCounts v//'"), emit: versions_subread, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def raw_args   = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def bam_inputs = bams instanceof List ? bams.join(' ') : "${bams}"

    // Normalize -t exon,CDS to single valid feature type (-t exon)
    def clean_args = raw_args.replaceAll(/-t\s+exon,CDS\b/, '-t exon')

    // Clean redundant -s flag from ext.args if present to prioritize meta.strandedness
    clean_args = clean_args.replaceAll(/(?:^|\s)-s\s*\d+(?:\s|$)/, ' ').trim()

    // Normalize paired-end / single-end flags
    def is_single_end = meta.single_end ? true : false
    if (is_single_end) {
        clean_args = clean_args.replaceAll(/(?:^|\s)-p(?:\s|$)/, ' ').trim()
    }
    def paired_end = (!is_single_end && !clean_args.contains('-p')) ? '-p' : ''
    clean_args = clean_args.replaceAll(/\s+/, ' ').trim()

    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    }
    else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }

    """
    featureCounts \\
        ${clean_args} \\
        ${paired_end} \\
        -T ${task.cpus} \\
        -a ${annotation} \\
        -s ${strandedness} \\
        -o ${prefix}.featureCounts.tsv \\
        ${bam_inputs}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.featureCounts.tsv
    touch ${prefix}.featureCounts.tsv.summary
    """
}
