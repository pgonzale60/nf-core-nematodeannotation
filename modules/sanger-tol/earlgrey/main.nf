process EARLGREY {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/earlgrey:6.1.1--h9948957_0' :
        'quay.io/biocontainers/earlgrey:6.1.1--h9948957_0' }"

    input:
    tuple val(meta), path(fasta)
    val rhabditida_clade

    output:
    // Paths verified against what EarlGrey 6.1.1 actually writes. The previous globs
    // pointed at the <prefix>_EarlGrey root, where EarlGrey writes NONE of these: the
    // merged calls and the RepeatLandscape live in <prefix>_mergedRepeats/ and
    // <prefix>_RepeatLandscape/, while the soft-masked genome and the summary
    // quantifications land in <prefix>_summaryFiles/ -- which EarlGrey itself announces
    // as its standard output location. That mismatch made emit:masked_fasta -- the only
    // non-optional emit -- unsatisfiable even on a fully successful run, so this task
    // could only ever fail, or be silently ignored under errorStrategy='ignore'.
    tuple val(meta), path("*_EarlGrey/*_summaryFiles/*.softmasked.fasta")             , emit: masked_fasta
    tuple val(meta), path("*_EarlGrey/*_summaryFiles/*.filteredRepeats.gff3")         , emit: gff3       , optional: true
    tuple val(meta), path("*_EarlGrey/*_summaryFiles/*.filteredRepeats.bed")          , emit: bed        , optional: true
    tuple val(meta), path("*_EarlGrey/*_summaryFiles/*.highLevelCount.txt")           , emit: summary    , optional: true
    tuple val(meta), path("*_EarlGrey/*_summaryFiles/*_combined_library.fasta")       , emit: te_library , optional: true
    tuple val(meta), path("*_EarlGrey/*_summaryFiles/*_divergence_summary_table.tsv") , emit: divergence , optional: true
    tuple val(meta), path("*_EarlGrey")                                              , emit: outdir
    tuple val("${task.process}"), val('earlgrey'), eval("earlGrey -h 2>&1 | grep -i 'earlgrey version' | head -1 | sed 's/.*version //; s/ //g'"), topic: versions, emit: versions_earlgrey

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '-e'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def clade  = rhabditida_clade ?: 'rhabditida'

    """
    # Run earlGrey non-interactively:
    # -d yes ensures soft-masked genome is generated at the end
    # Helitron search is controlled by args: the config passes -e yes (or -e no).
    #
    # --- PRE-FLIGHT -------------------------------------------------------------
    # EarlGrey is a monolith. In the 2026-10-01 run it spent ~24 hours completing
    # RepeatModeler, RepeatMasker, Strainer and library curation, then died in the
    # defragmentation stage because rmOutToGFF3.pl could not compile. That is a
    # five-second check. Everything here runs BEFORE the expensive work, so a
    # configuration defect costs seconds rather than a day.
    fail=0
    command -v earlGrey >/dev/null 2>&1 || { echo "FATAL: earlGrey not on PATH"; fail=1; }

    # rmOutToGFF3.pl ships without CrossmatchSearchEngine.pm on @INC in this
    # biocontainer AND in the equivalent conda env: the module is installed to
    # share/RepeatMasker, which perl never searches. Repair rather than rediscover it.
    for d in /opt/conda/share/RepeatMasker /usr/local/share/RepeatMasker /usr/share/RepeatMasker; do
        if [ -f "\$d/CrossmatchSearchEngine.pm" ]; then
            export PERL5LIB="\$d:\${PERL5LIB:-}"
        fi
    done

    printf '   SW   perc perc perc  query  position in query   matching  repeat  position in repeat\\n' > /tmp/earlgrey_preflight.out
    printf '    1    0.0  0.0  0.0  ctg_1  1  10  (0) +  Simple  1  10  (0)\\n' >> /tmp/earlgrey_preflight.out
    if ! rmOutToGFF3.pl /tmp/earlgrey_preflight.out > /dev/null 2> /tmp/earlgrey_preflight.err; then
        echo "FATAL: rmOutToGFF3.pl cannot run, so EarlGrey would fail at defragmentation."
        echo "       CrossmatchSearchEngine.pm is not on @INC (PERL5LIB=\${PERL5LIB:-unset})."
        head -3 /tmp/earlgrey_preflight.err
        fail=1
    fi

    case "${args}" in
        *"-e yes"*)
            command -v heliano >/dev/null 2>&1 || { echo "FATAL: -e yes requested but heliano is not in this container"; fail=1; }
            ;;
    esac

    [ "\$fail" -eq 0 ] || exit 1
    echo "preflight OK (PERL5LIB=\${PERL5LIB:-unset})"

    # Seed from prior run directory if configured (allows in-workflow resumption)
    if [ -n "${params.earlgrey_seed_dir ?: ''}" ] && [ -d "${params.earlgrey_seed_dir}" ]; then
        echo "Seeding EarlGrey working directory from: ${params.earlgrey_seed_dir}"
        target_dir="\$(basename "${params.earlgrey_seed_dir}")"
        if [ ! -d "\${target_dir}" ]; then
            cp -al "${params.earlgrey_seed_dir}" ./
        fi
    fi

    PERL5LIB="\${PERL5LIB:-}" \
    earlGrey \
        -g "${fasta}" \
        -s "${prefix}" \
        -o . \
        -r "${clade}" \
        -t ${task.cpus} \
        -d yes \
        ${args} </dev/null

    # --- POST-RUN CONTRACT ASSERTION -------------------------------------------
    # EarlGrey exits 0 before its deliverables necessarily exist, and the outputs
    # declared by this module previously pointed at paths it never writes -- so a
    # broken run could look successful. Assert the contract before returning.
    ok=0
    for f in "${prefix}_EarlGrey/${prefix}_summaryFiles/${prefix}.softmasked.fasta" \\
             "${prefix}_EarlGrey/${prefix}_summaryFiles/${prefix}.filteredRepeats.gff"; do
        if [ -s "\$f" ]; then echo "  contract OK: \$f"; else echo "  FATAL: missing or empty \$f"; ok=1; fi
    done
    if [ "\$ok" -ne 0 ]; then
        echo "FATAL: EarlGrey returned success without producing its declared deliverables."
        echo "       The intermediate tree is retained in this task directory; the staged"
        echo "       resume procedure is recorded in AGENT2_FINDINGS_AND_DECISIONS.md."
        exit 1
    fi

    # Normalize GFF extension if output as .gff.
    # EarlGrey writes <prefix>_EarlGrey/<prefix>_summaryFiles/<prefix>.filteredRepeats.gff
    # and nothing at the <prefix>_EarlGrey root, so the original condition here -- which
    # tested the root -- never fired, and the declared .gff3 output never existed.
    te_gff="${prefix}_EarlGrey/${prefix}_summaryFiles/${prefix}.filteredRepeats.gff"
    if [ -f "\${te_gff}" ] && [ ! -f "\${te_gff}3" ]; then
        cp "\${te_gff}" "\${te_gff}3"
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}_EarlGrey"
    touch "${prefix}_EarlGrey/${prefix}_filteredRepeats.softmasked.fasta"
    touch "${prefix}_EarlGrey/${prefix}_filteredRepeats.gff3"
    touch "${prefix}_EarlGrey/${prefix}_filteredRepeats.bed"
    touch "${prefix}_EarlGrey/${prefix}_summary.txt"
    """
}
