process TRF {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trf:4.09.1--h031d066_6' :
        'quay.io/biocontainers/trf:4.09.1--h031d066_6' }"

    input:
    tuple val(meta), path(fasta)
    val match
    val mismatch
    val delta
    val pm
    val pi
    val minscore
    val maxperiod

    output:
    tuple val(meta), path("*.dat"), emit: dat
    tuple val("${task.process}"), val('trf'), eval("strings \$(which trf) | grep -E '^[0-9]+(\\.[0-9]+)+' | head -n 1"), topic: versions, emit: versions_trf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: '-d -h'
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def m         = match ?: 2
    def mm        = mismatch ?: 7
    def d         = delta ?: 7
    def prob_m    = pm ?: 80
    def prob_i    = pi ?: 10
    def min_s     = minscore ?: 50
    def max_p     = maxperiod ?: 500

    """
    # TRF exits 20 when no repeats found in a window; || true lets the
    # script continue so the .dat normalisation below can always run.
    trf \\
        "${fasta}" \\
        ${m} \\
        ${mm} \\
        ${d} \\
        ${prob_m} \\
        ${prob_i} \\
        ${min_s} \\
        ${max_p} \\
        ${args} || true

    # TRF outputs <fasta>.<m>.<mm>.<d>.<prob_m>.<prob_i>.<min_s>.<max_p>.dat
    # Normalize output name to ${prefix}.dat
    DAT_FILE=\$(ls *.dat 2>/dev/null | head -n 1 || true)
    if [ -n "\$DAT_FILE" ] && [ "\$DAT_FILE" != "${prefix}.dat" ]; then
        mv "\$DAT_FILE" "${prefix}.dat"
    elif [ -z "\$DAT_FILE" ]; then
        touch "${prefix}.dat"
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dat
    """
}
