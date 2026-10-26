process INTERPROSCAN {
    tag "${meta.id}"
    label 'process_high'

    // InterProScan requires extensive data files; bind or reference host/shared installation
    input:
    tuple val(meta), path(fasta)
    path ipr_dir

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.gff3"), emit: gff3, optional: true
    tuple val(meta), path("${prefix}.xml"), emit: xml, optional: true
    tuple val("${task.process}"), val('interproscan'), val('5.63-95.0'), emit: versions_interproscan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-f tsv,gff3,xml -iprlookup -goterms -pathways -pa -dp'
    prefix   = task.ext.prefix ?: "${meta.id}.iprscan"
    """
    # 1. Sanitize input: Strip '*' and convert internal/frameshift '.' stop characters to 'X' (IUPAC valid)
    sed '/^[^>]/ { s/\\*//g; s/\\./X/g; }' "${fasta}" > "${meta.id}.clean.fasta"

    # 2. Run InterProScan using container or local/shared installation
    IPR_EXEC="interproscan.sh"
    if [ -x "${ipr_dir}/interproscan.sh" ]; then
        IPR_EXEC="${ipr_dir}/interproscan.sh"
    elif [ -n "\${JAVA_HOME:-}" ]; then
        export PATH="\${JAVA_HOME}/bin:\${PATH}"
    fi

    mkdir -p tmp results

    # When multiple formats are requested (-f tsv,gff3,xml), InterProScan requires -b (output base name) instead of -o
    "\${IPR_EXEC}" \\
        -i "${meta.id}.clean.fasta" \\
        ${args} \\
        -cpu ${task.cpus} \\
        -T tmp \\
        -b "${prefix}"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.iprscan"
    """
    touch "${prefix}.tsv"
    touch "${prefix}.gff3"
    touch "${prefix}.xml"
    """
}
