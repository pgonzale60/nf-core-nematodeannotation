process VALIDATE_RNASEQ_ALIGNMENT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c2f02ceb8cbe5aa37930dcbc0de6f1d07dfb292178c7e780c44ab2f556bf829f/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--8a3ddf262165ba23'}"

    input:
    tuple val(meta), path(alignment), path(index)
    tuple val(meta_fai), path(fai)

    output:
    tuple val(meta), path(alignment), path(index), env('IS_HISAT2'), env('IS_COORDINATE_VALID'), emit: evaluation
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # 1. Structural integrity check
    if ! samtools quickcheck -q "${alignment}"; then
        echo "ERROR: Alignment file '${alignment}' failed samtools quickcheck (corrupted or truncated)." >&2
        exit 1
    fi

    # 2. Alignment provenance check (HISAT2)
    if samtools view -H "${alignment}" | grep -E '^@PG.*(ID:hisat2|PN:hisat2|hisat2-align)' >/dev/null 2>&1; then
        echo "Found HISAT2 provenance in header of ${alignment}"
        IS_HISAT2="true"
    else
        echo "HISAT2 provenance NOT found in header of ${alignment}"
        IS_HISAT2="false"
    fi

    # 3. Reference sequence & coordinate validation against assembly .fai
    samtools view -H "${alignment}" | awk '/^@SQ/ {for(i=1;i<=NF;i++){if(\$i~/^SN:/){sn=substr(\$i,4)}if(\$i~/^LN:/){ln=substr(\$i,4)}} print sn"\\t"ln}' | sort > bam_sq.txt
    awk '{print \$1"\\t"\$2}' "${fai}" | sort > fai_sq.txt

    # Verify that reference sequences in BAM header exist in target assembly FAI with identical lengths
    mismatches=\$(comm -23 bam_sq.txt fai_sq.txt | wc -l | tr -d ' ')
    if [ "\${mismatches}" -gt 0 ]; then
        echo "ERROR: Reference sequence mismatch detected in '${alignment}'." >&2
        echo "The following contigs/lengths in the BAM header do not match the target assembly:" >&2
        comm -23 bam_sq.txt fai_sq.txt >&2
        IS_COORDINATE_VALID="false"
        exit 1
    else
        echo "Reference sequence dictionary in '${alignment}' matches target assembly."
        IS_COORDINATE_VALID="true"
    fi
    """

    stub:
    """
    IS_HISAT2="true"
    IS_COORDINATE_VALID="true"
    """
}
