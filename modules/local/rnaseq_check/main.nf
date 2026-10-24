process CHECK_HISAT2_ALIGNMENT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c2f02ceb8cbe5aa37930dcbc0de6f1d07dfb292178c7e780c44ab2f556bf829f/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--8a3ddf262165ba23'}"

    input:
    tuple val(meta), path(alignment), path(index)
    path fasta

    output:
    tuple val(meta), path(alignment), env('IS_HISAT2'), emit: evaluation
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def ref_arg = fasta ? "--reference ${fasta}" : ""
    """
    if samtools view -H ${ref_arg} "${alignment}" | grep -E '^@PG.*(ID:hisat2|PN:hisat2|hisat2-align)' >/dev/null 2>&1; then
        echo "Found HISAT2 provenance in header of ${alignment}"
        IS_HISAT2="true"
    else
        echo "HISAT2 provenance NOT found in header of ${alignment}"
        IS_HISAT2="false"
    fi
    """

    stub:
    """
    IS_HISAT2="true"
    """
}
