process SELECTBESTMODEL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/97/97c07dca6975c8aa84249314cc7aedb72cf651b77f4941a3cf98104e84eb8b28/data'
        : 'community.wave.seqera.io/library/busco_numpy:3ecb8547e6c97dd1'}"

    input:
    tuple val(meta), val(labels), path(gtfs), path(aas)
    val val_busco_lineage
    path busco_lineages_path

    output:
    tuple val(meta), path("${prefix}.gtf"),          emit: representative_gtf
    tuple val(meta), path("${prefix}.aa"),           emit: representative_aa
    tuple val(meta), path("${prefix}_decision.tsv"), emit: decision_tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.representative"
    def lineage_dir_arg = busco_lineages_path ? "--lineage-dir ${busco_lineages_path}" : ""
    def label_list = labels instanceof List ? labels : [labels]
    def gtf_list   = gtfs instanceof List ? gtfs : [gtfs]
    def aa_list    = aas instanceof List ? aas : [aas]
    def cand_args  = []
    [label_list, gtf_list, aa_list].transpose().each { l, g, a ->
        cand_args.add("--candidate ${l} ${g} ${a}")
    }
    def cand_str = cand_args.join(' ')
    """
    # Tool: select_best_annotation.py | Script-SHA: 71c6fa2d5c9b9934
    # Fix Augustus environment variables for Apptainer if present
    ENV_AUGUSTUS=/opt/conda/etc/conda/activate.d/augustus.sh
    if [ -f "\$ENV_AUGUSTUS" ]; then
        source "\$ENV_AUGUSTUS"
    fi

    python3 ${moduleDir}/resources/usr/bin/select_best_annotation.py \\
        --lineage "${val_busco_lineage}" \\
        ${lineage_dir_arg} \\
        --cpus ${task.cpus} \\
        --prefix "${prefix}" \\
        --outdir . \\
        ${cand_str}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.representative"
    """
    touch "${prefix}.gtf"
    touch "${prefix}.aa"
    touch "${prefix}_decision.tsv"
    """
}
