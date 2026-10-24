process BLASTTOBED {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(blast_tsv)

    output:
    tuple val(meta), path("*.sorted.bed"), emit: bed
    tuple val("${task.process}"), val('awk'), eval("awk 2>&1 | head -n 1 | sed -e 's/.*BusyBox v//; s/ .*//; s/GNU Awk //'"), topic: versions, emit: versions_awk

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_gz = blast_tsv.getExtension() == "gz"
    def cat_cmd = is_gz ? "gzip -cd" : "cat"

    """
    ${cat_cmd} "${blast_tsv}" | awk -F'\\t' 'BEGIN {OFS="\\t"} {
        if (\$0 !~ /^#/) {
            sstart = \$9 + 0;
            send   = \$10 + 0;
            start  = (sstart < send) ? sstart - 1 : send - 1;
            end    = (sstart < send) ? send : sstart;
            strand = (sstart <= send) ? "+" : "-";
            if (end > start && start >= 0) {
                print \$2, start, end, \$1, int(\$3), strand;
            }
        }
    }' | sort -k1,1 -k2,2n > "${prefix}.sorted.bed"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.sorted.bed"
    """
}
