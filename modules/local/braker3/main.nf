process BRAKER3 {
    tag "${meta.id}"
    label 'process_high'

    // Re. Conda from the BRAKER team:
    // Warning: installing GeneMark-ETP for BRAKER in conda environments has lead to multiple problems reported by users (Issues!).
    // We can not offer support for conda installations. Please use the singularity image instead.
    container "docker.io/teambraker/braker3:v3.0.7.5"

    input:
    tuple val(meta), path(fasta)
    path bam
    path rnaseq_sets_dirs
    path rnaseq_sets_ids
    path proteins
    path hintsfile

    output:
    tuple val(meta), path("$prefix/braker.gtf")                     , emit: gtf
    tuple val(meta), path("$prefix/braker.codingseq")               , emit: cds
    tuple val(meta), path("$prefix/braker.aa")                      , emit: aa
    tuple val(meta), path("$prefix/Augustus/augustus.hints.gtf")     , emit: augustus_gtf, optional: true
    tuple val(meta), path("$prefix/Augustus/augustus.hints.aa")      , emit: augustus_aa , optional: true
    tuple val(meta), path("$prefix/GeneMark-*/genemark.gtf")         , emit: genemark_gtf, optional: true
    tuple val(meta), path("$prefix/braker.log")                     , emit: log
    tuple val(meta), path("$prefix/hintsfile.gff")                  , emit: hintsfile   , optional: true
    tuple val(meta), path("$prefix/braker.gff3")                    , emit: gff3        , optional: true
    tuple val(meta), path("$prefix/what-to-cite.txt")               , emit: citations
    tuple val("${task.process}"), val('braker3'), eval("braker.pl --version 2>/dev/null | sed 's/braker.pl version //'"), emit: versions_braker3, topic: versions
    tuple val("${task.process}"), val('augustus'), eval("augustus --version |& sed -n 's/AUGUSTUS (\\(.*\\)) is a gene .*/\\1/p'"), emit: versions_augustus, topic: versions
    tuple val("${task.process}"), val('genemark-etp'), eval("gmetp.pl |& sed -n 's/ETP version \\(.*\\)/\\1/p'"), emit: versions_genemarketp, topic: versions
    tuple val("${task.process}"), val('prothint'), eval("prothint.py --version |& sed 's/prothint.py //1'"), emit: versions_prothint, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args             ?: ''
    prefix          = task.ext.prefix           ?: "${meta.id}"
    def rna_ids     = rnaseq_sets_ids           ? "--rnaseq_sets_ids=$rnaseq_sets_ids"      : ''
    def rna_dirs    = rnaseq_sets_dirs          ? "--rnaseq_sets_dirs=$rnaseq_sets_dirs"    : ''
    def bam_arg     = bam                       ? "--bam=$bam"                              : ''
    def prot_arg    = proteins                  ? "--prot_seq=$proteins"                    : ''
    def hints       = hintsfile                 ? "--hints=$hintsfile"                      : ''
    def new_species = args.contains('--species')? ''                                        : '--species new_species'
    """
    # --- PRE-FLIGHT -------------------------------------------------------------
    # BRAKER3 is a monolith that runs for hours to days. The failures that actually
    # waste that time are configuration ones -- a BAM without an index, an empty
    # protein set, an absent GeneMark key -- and every one of them is detectable in
    # seconds. Detect them here, before the run rather than deep inside it.
    fail=0
    for t in braker.pl augustus gmetp.pl prothint.py samtools; do
        command -v \$t >/dev/null 2>&1 || { echo "FATAL: \$t not on PATH"; fail=1; }
    done

    if [ -n "${bam_arg}" ]; then
        samtools quickcheck -v "$bam" || { echo "FATAL: BAM fails samtools quickcheck: $bam"; fail=1; }
        samtools idxstats "$bam" > /dev/null 2>&1 || { echo "FATAL: BAM has no usable index (.bai missing?): $bam"; fail=1; }
    fi

    if [ -n "${prot_arg}" ]; then
        n_prot=\$(grep -c '^>' "$proteins" || true)
        [ "\$n_prot" -gt 0 ] || { echo "FATAL: protein FASTA has no records: $proteins"; fail=1; }
        echo "  proteins : \$n_prot records"
    fi

    n_genome=\$(grep -c '^>' $fasta || true)
    [ "\$n_genome" -gt 0 ] || { echo "FATAL: genome FASTA has no records: $fasta"; fail=1; }
    echo "  genome   : \$n_genome sequences"

    if [ ! -e "\$HOME/.gm_key" ] && [ ! -e /root/.gm_key ] && [ ! -e /opt/GeneMark-ETP/gm_key ]; then
        echo "WARNING: no GeneMark key found (~/.gm_key, /root/.gm_key, /opt/GeneMark-ETP/gm_key)."
        echo "         gmetp.pl will fail if this invocation requires one."
    fi

    [ "\$fail" -eq 0 ] || exit 1
    echo "preflight OK"

    AUG_CONF="\${AUGUSTUS_CONFIG_PATH:-}"
    if [ -z "\$AUG_CONF" ] || [ ! -d "\$AUG_CONF" ]; then
        if [ -d "/opt/Augustus/config" ]; then
            AUG_CONF="/opt/Augustus/config"
        elif [ -d "/opt/conda/config" ]; then
            AUG_CONF="/opt/conda/config"
        else
            AUG_CONF="\$(dirname \$(which augustus))/../config"
        fi
    fi

    export PYTHONNOUSERSITE=1
    export TMPDIR="\${PWD}/tmp"
    mkdir -p "\$TMPDIR"

    cp -r "\$AUG_CONF" augustus_config
    chmod -R a+w augustus_config

    perl -p -e 's/^(>\\S+).*\$/\$1/' \\
        $fasta \\
        > ${prefix}.name.only.genome.masked.fasta

    braker.pl \\
        --genome ${prefix}.name.only.genome.masked.fasta \\
        $new_species \\
        --workingdir $prefix \\
        --AUGUSTUS_CONFIG_PATH "\$(pwd)/augustus_config" \\
        --threads $task.cpus \\
        $rna_ids \\
        $rna_dirs \\
        $bam_arg \\
        $prot_arg \\
        $hints \\
        $args

    # --- POST-RUN CONTRACT ASSERTION -------------------------------------------
    # braker.gtf is a declared, non-optional output, so Nextflow would catch its
    # absence -- but only once the whole run is over. Assert the useful outputs here
    # so a hollow success fails immediately and legibly instead.
    for f in "$prefix/braker.gtf" "$prefix/braker.aa" "$prefix/braker.log"; do
        [ -s "\$f" ] || { echo "FATAL: BRAKER3 produced no non-empty \$f"; exit 1; }
    done
    echo "contract OK: braker.gtf, braker.aa and braker.log are non-empty"
    """

    stub:
    def args        = task.ext.args                         ?: ''
    prefix          = task.ext.prefix                       ?: "${meta.id}"
    def rna_ids     = rnaseq_sets_ids                       ? "--rnaseq_sets_ids=$rnaseq_sets_ids"      : ''
    def hints       = hintsfile                             ? "--hints=$hintsfile"                      : ''
    def touch_hints = (rna_ids || bam || proteins || hints) ? "touch $prefix/hintsfile.gff"             : ''
    def touch_gff   = args.contains('--gff3')               ? "touch $prefix/braker.gff3"               : ''
    """
    mkdir "$prefix"

    touch "$prefix/braker.gtf"
    touch "$prefix/braker.codingseq"
    touch "$prefix/braker.aa"
    $touch_hints
    touch "$prefix/braker.log"
    touch "$prefix/what-to-cite.txt"
    $touch_gff
    """
}
